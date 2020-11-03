from datetime import datetime, timedelta
import pathlib
import shutil
import urllib.request
import urllib.error

from astropy.coordinates import SkyCoord, Longitude
from astropy.time import Time
import astropy.units as u
from astropy.wcs import WCS
from reproject import reproject_interp

import numpy as np
import matplotlib.pyplot as plt

from sunpy.net import vso
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map, make_fitswcs_header
from sunpy.coordinates import get_earth
import sunpy.sun.constants

from aiapy.calibrate import (update_pointing, fix_observer_location,
                             correct_degradation, normalize_exposure)
from aiapy.calibrate.util import get_correction_table

from time_helpers import start_of_day


# The directory in which maps are downloaded to and read from
map_dir = pathlib.Path('/Volumes/Work/Data/aia_new')
if not map_dir.exists():
    raise RuntimeError(f'Map directory {map_dir} does not exist')
correction_table = get_correction_table()


def map_path(dtime, wlen):
    """
    Get the path of an AIA map at a given date and wavelength.
    """
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_{wlen}_{datestr}.fits'


def synoptic_map_path(dtime, wlen):
    """
    Get the path of an AIA synoptic map at a given date and wavelength.
    """
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_{wlen}_synoptic_{datestr}.fits'


def download_start_of_day_map(dtime, wlen):
    """
    Download the first map available on a given date, at a given wavelength.
    """
    dtime = start_of_day(dtime)
    print(f'Fetching map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(wlen * u.Angstrom))
    result = Fido.search(*query)
    try:
        download_path = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No map available for {dtime}')
    download_path = pathlib.Path(download_path)
    shutil.move(download_path, map_path(dtime, wlen))


def load_start_of_day_map(dtime, wlen):
    """
    Load the first map available on a given date, at a given wavelength.
    """
    dtime = start_of_day(dtime)
    mappath = map_path(dtime, wlen)
    if not mappath.exists():
        download_start_of_day_map(dtime, wlen)

    print(f'Loading AIA {wlen} map for {dtime}')
    return Map(str(mappath))


def prep(m):
    """
    Prep an AIA map
    """
    print('Prepping map')
    if m.exposure_time <= 0 * u.s:
        raise RuntimeError('Exposure time <= 0')
    m = update_pointing(m)
    m = fix_observer_location(m)
    m = correct_degradation(m, correction_table=correction_table)
    m = normalize_exposure(m)
    return m


def synop_header(shape_out, dtime):
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime,
                         observer='earth')
    header = make_fitswcs_header(
        shape_out, frame_out,
        scale=[180 / shape_out[0],
               360 / shape_out[1]] * u.deg / u.pix,
        projection_code="CAR")
    return header


def helioproj_header(shape_out, dtime):
    sun_width_arscec = 2000
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame='helioprojective',
                         obstime=dtime, observer='earth')
    header = make_fitswcs_header(
        shape_out, frame_out,
        scale=[sun_width_arscec / shape_out[0],
               sun_width_arscec / shape_out[1]] * u.arcsec / u.pix)
    return header


def synop_reproject(m, shape_out, wlen):
    synop_map_path = synoptic_map_path(m.date, wlen)
    if not synop_map_path.exists():
        print(f'Reprojecting AIA {wlen} map')
        m = prep(m)
        m.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
        header = synop_header(shape_out, m.date)
        wcs = WCS(header)
        wcs.heliographic_observer = m.observer_coordinate
        with np.errstate(invalid='ignore'):
            array, footprint = reproject_interp(m, wcs, shape_out=shape_out)
        new_map = Map((array, header))
        new_map.save(str(synop_map_path))

    new_map = Map(synop_map_path)
    new_map.plot_settings = m.plot_settings
    return new_map


def long_weights(longs, l0):
    """
    Weights to use when adding synoptic maps.

    Parameters
    ----------
    longs :
        The longitude coordinates of each pixel in a map.
    l0 :
        The observer longitude.
    """
    dcenterlong = (longs - l0 + 180) % 360 - 180
    weights = np.exp(-(dcenterlong / 15)**2)
    weights[weights < 0] = 0
    return weights / np.nanmax(weights)


def create_synoptic_map(endtime, wlen):
    """
    Create a synoptic map, using 27 daily SDO/AIA maps ending on the
    endtime given. Note that the maps are taken from the start of each day.

    Returns
    -------
    sunpy.map.Map : synoptic map
    """
    shape = [720, 1440]
    data = np.zeros(shape)
    weight_sum = np.zeros(shape)
    nmaps = 27
    recent_time = None
    for i in range(nmaps):
        dtime = endtime - timedelta(days=i)
        try:
            aia_map = load_start_of_day_map(dtime, wlen)
            aia_synop_map = synop_reproject(aia_map, shape, wlen)
        except RuntimeError as e:
            print(e)
            continue

        if recent_time is None:
            recent_time = dtime.strftime('%Y-%m-%dT%H:%M:%S')

        # Create weights
        coord = sunpy.map.all_coordinates_from_map(aia_synop_map)
        longs = coord.lon.to(u.deg).value
        l0 = sunpy.coordinates.sun.L0(dtime).to(u.deg).value
        weights = long_weights(longs, l0)

        aia_data = aia_synop_map.data
        aia_data[np.isnan(aia_data)] = 0
        data += (aia_data * weights)
        weight_sum += weights

    weight_sum[weight_sum == 0] = np.nan
    data /= weight_sum

    meta = aia_synop_map.meta
    meta['date-obs'] = recent_time

    synop_map = Map((data, meta))
    synop_map.plot_settings = aia_synop_map.plot_settings
    return synop_map


def aia_fov(dtime):
    l0 = sunpy.coordinates.sun.L0(dtime)
    bounds = Longitude([l0 - 90 * u.deg, l0 + 90 * u.deg])
    return bounds


if __name__ == '__main__':
    map = create_synoptic_map(datetime(2018, 11, 10), 193)
    # Norm the data
    # data = map.data
    # data = map.plot_settings['norm'](data)

    # map = Map((data, map.meta))
    datestr = datetime.now().strftime('%Y%m%d')
    map.save(f'aia193_synoptic_latest_{datestr}.fits')

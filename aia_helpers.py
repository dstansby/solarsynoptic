from datetime import datetime, timedelta
import pathlib
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

from time_helpers import start_of_day


map_dir = pathlib.Path('/Users/dstansby/Data/aia')
map_dir.mkdir(exist_ok=True, parents=True)


def map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_193_{datestr}.fits'


def synoptic_map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_193_synoptic_{datestr}.fits'


def download_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    print(f'Fetching map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(193 * u.Angstrom))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No map available for {dtime}')
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def load_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    print(f'Loading AIA map for {dtime}')
    return Map(str(mappath))


def synop_header(shape_out, dtime):
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime)
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


def synop_reproject(m, shape_out):
    synop_map_path = synoptic_map_path(m.date.to_datetime())
    if not synop_map_path.exists():
        m.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
        header = synop_header(shape_out, m.date)
        array, footprint = reproject_interp(m, WCS(header),
                                            shape_out=shape_out)
        new_map = Map((array, header))
        new_map.save(str(synop_map_path))

    print(f'Loading {synop_map_path}')
    new_map = Map(synop_map_path)
    new_map.plot_settings = m.plot_settings
    return new_map


def create_synoptic_map(endtime):
    """
    Create an AIA synoptic map, using 27 daily AIA 193 maps ending on the
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
        aia_map = load_start_of_day_map(dtime)
        aia_synop_map = synop_reproject(aia_map, shape)

        if recent_time is None:
            recent_time = dtime.strftime('%Y-%m-%dT%H:%M:%S')

        # Create weights
        coord = sunpy.map.all_coordinates_from_map(aia_synop_map)
        longs = coord.lon.to(u.deg).value
        l0 = sunpy.coordinates.sun.L0(dtime).to(u.deg).value
        dcenterlong = (longs - l0 + 180) % 360 - 180
        weights = np.exp(-(dcenterlong / 10)**2)
        weights[weights < 0] = 0

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


def stonyhurst_reproject(synoptic_map, dtime):
    '''
    Reproject a synoptic map into an Earth facing view
    (ie. a heliographic Stonyhurst frame).
    '''
    shape_out = [1024, 1024]
    header = helioproj_header(shape_out, dtime)
    wcs = WCS(header)
    # Now have to manually set the observer coordinate
    wcs.heliographic_observer = get_earth(dtime)
    # array, footprint = reproject_interp(synoptic_map, wcs,
    #                                     shape_out=shape_out)
    array = np.random.rand(shape_out[0], shape_out[1])
    new_map = Map((array, header))
    new_map.plot_settings = synoptic_map.plot_settings
    return new_map


def aia_fov(dtime):
    l0 = sunpy.coordinates.sun.L0(dtime)
    bounds = Longitude([l0 - 90 * u.deg, l0 + 90 * u.deg])
    return bounds


if __name__ == '__main__':
    map = create_synoptic_map(datetime.now())
    # Norm the data
    data = map.data
    data = map.plot_settings['norm'](data)

    map = Map((data, map.meta))
    datestr = datetime.now().strftime('%Y%m%d')
    map.save(f'aia193_synoptic_latest_{datestr}.fits')

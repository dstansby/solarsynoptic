import pathlib
import shutil
from datetime import datetime, timedelta

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.sun.constants
from aiapy.calibrate import (
    correct_degradation,
    normalize_exposure,
    update_pointing,
)
from aiapy.calibrate.util import get_correction_table
from astropy.coordinates import Longitude, SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp
from sunpy.map import Map, make_fitswcs_header
from sunpy.net import Fido
from sunpy.net import attrs as a

from time_helpers import start_of_day

# The directory in which maps are downloaded to and read from
map_dir = pathlib.Path('/Volumes/Work/Data/aia')
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


def download_nrt(dtime, wlen):
    import parfive
    dl = parfive.Downloader(max_conn=1)
    url = (f"http://jsoc2.stanford.edu/data/aia/synoptic/nrt/"
           f"{dtime.year}/{dtime.month:02}/{dtime.day:02}/"
           f"H0000/AIA{dtime.year}{dtime.month:02}{dtime.day:02}_"
           f"000000_0{wlen}.fits")
    dl.enqueue_file(url, filename=map_path(dtime, wlen))
    res = dl.download()
    if len(res.errors):
        print(res.errors)
        raise RuntimeError('Download failed')
    download_path = pathlib.Path(res[0])
    shutil.move(download_path, map_path(dtime, wlen))


def download_start_of_day_map(dtime, wlen):
    """
    Download the first map available on a given date, at a given wavelength.
    """
    if dtime > datetime.now():
        raise RuntimeError(f'No map available for {dtime}')
    elif dtime < datetime.now() - timedelta(days=7):
        dtime = start_of_day(dtime)
        print(f'Fetching map for {dtime}')
        query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
                 a.Instrument('AIA'),
                 a.Wavelength(wlen * u.Angstrom))
        result = Fido.search(*query)
        try:
            download_path = Fido.fetch(result[0, 0])[0]
        except IndexError:
            raise ValueError(f'Failed to donwload map for {dtime}')
        download_path = pathlib.Path(download_path)
        shutil.move(download_path, map_path(dtime, wlen))
    else:
        download_nrt(dtime, wlen)


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
    Prep an AIA map.

    This runs (in order):
        - `aiapy.calibrate.update_pointing`
        - `aiapy.calibrate.fix_observer_location`
        - `aiapy.caibrate.correct_degradation`
        - `aiapy.calibrate.normalize_exposure`
    """
    print('Prepping map')
    if m.exposure_time <= 0 * u.s:
        raise RuntimeError('Exposure time <= 0')
    m = update_pointing(m)
    # m = fix_observer_location(m)
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
    nmaps = 28
    recent_time = None
    for i in range(nmaps):
        dtime = endtime - timedelta(days=i)
        try:
            aia_map = load_start_of_day_map(dtime, wlen)
            aia_synop_map = synop_reproject(aia_map, shape, wlen)
        except (RuntimeError, KeyError, ValueError) as e:
            print('\U0001F6A8 ' + str(e))
            continue

        if recent_time is None:
            recent_time = dtime.strftime('%Y-%m-%dT%H:%M:%S')

        # Create weights
        coord = sunpy.map.all_coordinates_from_map(aia_synop_map)
        longs = coord.lon.to(u.deg).value
        l0 = sunpy.coordinates.sun.L0(dtime).to(u.deg).value
        weights = long_weights(longs, l0)

        aia_data = aia_synop_map.data
        # Cast missing data to zero for now, to avoid adding NaNs to the total
        aia_data[~np.isfinite(aia_data)] = 0

        data += (aia_data * weights)
        weight_sum += weights

    data /= weight_sum
    data[data == 0] = np.nan

    meta = aia_synop_map.meta
    meta['date-obs'] = recent_time
    meta['instrume'] = 'AIA'    # Set so sunpy recognises this as an AIA map
    meta['telescop'] = 'SDO'
    meta['wavelnth'] = str(wlen)
    meta['waveunit'] = 'Angstrom'

    synop_map = Map((data, meta))
    synop_map.plot_settings = aia_synop_map.plot_settings
    return synop_map


def aia_fov(dtime):
    l0 = sunpy.coordinates.sun.L0(dtime)
    bounds = Longitude([l0 - 90 * u.deg, l0 + 90 * u.deg])
    return bounds


if __name__ == '__main__':
    end_date = datetime.now()
    map = create_synoptic_map(end_date, 193)
    # Norm the data
    # data = map.data
    # data = map.plot_settings['norm'](data)

    # map = Map((data, map.meta))
    datestr = end_date.strftime('%Y%m%d')
    map.save(f'aia193_synoptic_latest_{datestr}.fits')

    map.plot(vmin=0, vmax=2000)
    plt.gcf().savefig(f'aia193_synoptic_latest_png_{datestr}.png')

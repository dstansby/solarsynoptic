from datetime import datetime, timedelta
import pathlib

import astropy.units as u
from sunpy.map import Map
from sunpy.net import attrs as a
from sunpy.time import parse_time

from . import helpers
# TODO: change this to get the sunpy directory from sunpy
map_dir = pathlib.Path('/Users/dstansby/sunpy/solarsynoptic/raw')

__all__ = ['aia_start_of_day_map']


def aia_start_of_day_map(dtime, wlen):
    """
    Download the first AIA map available on a given date, at a given wavelength.

    Parameters
    ----------
    dtime
    wlen : astropy.units.Quantity
        Wavelength of interest.
    """
    dtime = parse_time(dtime).to_datetime()
    dtime = helpers.start_of_day(dtime)

    # Download from JSOC if older than 14 days; otherwise directly
    # get an NRT map
    if map_path(dtime, wlen).exists():
        pass
    elif (dtime < datetime.now() - timedelta(days=13)):
        query = a.Instrument('AIA'), a.Wavelength(wlen)
        mappath = helpers.start_of_day_map(dtime, *query)
        mappath = pathlib.Path(mappath)
        mappath.replace(map_path(dtime, wlen))
    else:
        import parfive
        dl = parfive.Downloader(max_conn=1)
        wlen_int = int(wlen.to_value(u.Angstrom))
        url = (f"http://jsoc2.stanford.edu/data/aia/synoptic/nrt/"
               f"{dtime.year}/{dtime.month:02}/{dtime.day:02}/"
               f"H0000/AIA{dtime.year}{dtime.month:02}{dtime.day:02}_"
               f"000000_0{wlen_int}.fits")
        dl.enqueue_file(url, filename=map_path(dtime, wlen))
        res = dl.download()
        if len(res.errors):
            print(res.errors)
            raise RuntimeError('Download failed')
        mappath = pathlib.Path(res[0])
        mappath.replace(map_path(dtime, wlen))

    return Map(map_path(dtime, wlen))


def map_path(dtime, wlen):
    datestr = dtime.strftime('%Y%m%d')
    wlen = int(wlen.to_value(u.Angstrom))
    return map_dir / f'aia_{wlen}_nrt_{datestr}.fits'

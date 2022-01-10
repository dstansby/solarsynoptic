import pathlib
from datetime import timedelta

import astropy.units as u
from parfive import Downloader
from sunpy.map import Map
from sunpy.time import TimeRange

from . import helpers

# TODO: change this to get the sunpy directory from sunpy
map_dir = pathlib.Path('/Users/dstansby/sunpy/solarsynoptic/raw')


__all__ = ['stereo_start_of_day_map']


def map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'euvi_195_{datestr}.fits'


def stereo_beacon_start_of_day(dtime):
    from sunpy.util.scraper import Scraper
    datestr = dtime.strftime('%Y%m%d')
    scrp = Scraper('https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/secchi/img'
                   f'/euvi/{datestr}/%Y%m%d_%H%M%S_n7euA.fts')
    urls = scrp.filelist(TimeRange(dtime, dtime + timedelta(days=1)))
    urls.sort()
    if not len(urls):
        raise RuntimeError(f'No EUVI beacon map available for {dtime}')
    dl = Downloader()
    have_195 = False
    for url in urls:
        dl.enqueue_file(url, path=map_path(dtime).parent)
        files = dl.download()
        if len(files.errors):
            raise RuntimeError(f'EUVI download failed for {dtime}')
        if Map(files[0]).wavelength != 195 * u.Angstrom:
            continue
        else:
            pathlib.Path(files[0]).replace(map_path(dtime))
            have_195 = True
            break
    if not have_195:
        raise RuntimeError('No 195 angstrom images')
    return map_path(dtime)


def download_start_of_day_map(dtime):
    dtime = helpers.start_of_day(dtime)
    '''query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('EUVI'))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:'''
    mappath = stereo_beacon_start_of_day(dtime)
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def stereo_start_of_day_map(dtime):
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    return Map(mappath)

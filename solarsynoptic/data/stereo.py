import pathlib

from parfive import Downloader
from sunpy.map import Map

from . import helpers
# TODO: change this to get the sunpy directory from sunpy
map_dir = pathlib.Path('/Users/dstansby/sunpy/data')


__all__ = ['stereo_start_of_day_map']


def map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'euvi_195_{datestr}.fits'


def stereo_beacon_start_of_day(dtime):
    datestr = dtime.strftime('%Y%m%d')
    for time in ['000530', '001615', '005530', '010530', '011530', '043530', ]:
        url = ('https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/secchi/img'
               f'/euvi/{datestr}/{datestr}_{time}_n7euA.fts')
        dl = Downloader()
        dl.enqueue_file(url, path=map_path(dtime).parent)
        files = dl.download()
        if len(files.errors):
            continue
        pathlib.Path(files[0]).replace(map_path(dtime))
        return
    raise RuntimeError(f'No EUVI beacon map available for {dtime}')


def download_start_of_day_map(dtime):
    dtime = helpers.start_of_day(dtime)
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('EUVI'))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        download_beacon(dtime)
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def stereo_start_of_day_map(dtime):
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    return Map(mappath)

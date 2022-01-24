from datetime import timedelta
from pathlib import Path

import astropy.units as u
from parfive import Downloader
from sunpy import config as sunpy_config
from sunpy.map import Map
from sunpy.net import attrs as a
from sunpy.time import TimeRange

from . import helpers

__all__ = ['stereo_start_of_day_map']


def _stereo_map_path(dtime, directory, wlen, beacon=''):
    datestr = dtime.strftime('%Y%m%d')
    wlen = int(wlen.to_value(u.Angstrom))
    return Path(directory) / f'euvi_{wlen}_{datestr}{beacon}.fits'


def stereo_beacon_start_of_day(dtime, directory, wlen=195*u.Angstrom):
    """
    Download the STEREO beacon file available on date *dtime* to *directory*.

    Returns
    -------
    pathlib.Path
        Path to downloaded file.

    Notes
    -----
    STEREO beacon files do not contain their wavelength in the filename,
    so this funciton downloads files one at a time and opens them to check
    if they have the right wavelength.
    """
    from sunpy.net import Scraper
    datestr = dtime.strftime('%Y%m%d')
    scrp = Scraper('https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/secchi/img'
                   f'/euvi/{datestr}/%Y%m%d_%H%M%S_n7euA.fts')
    urls = scrp.filelist(TimeRange(dtime, dtime + timedelta(days=1)))
    urls.sort()
    if not len(urls):
        raise RuntimeError(f'No EUVI beacon maps available on date {dtime}')

    dl = Downloader()
    have_right_wlen = False
    for url in urls:
        dl.enqueue_file(url, path=directory)
        files = dl.download()
        if len(files.errors):
            raise RuntimeError(f'EUVI download failed for {dtime}')
        if Map(files[0]).wavelength != wlen:
            continue
        else:
            have_right_wlen = True
            break
    if not have_right_wlen:
        raise RuntimeError(f'No {wlen} images found')
    return Path(files[0])


def stereo_start_of_day_map(dtime, wlen=195 * u.Angstrom, dl_path=None):
    """
    Download and load the first STEREO EUVI map available on a given date at a
    given wavelength.

    By default files are downloaded to the sunpy download directory.

    Parameters
    ----------
    dtime : astropy.time.Time
    wlen : astropy.units.Quantity
        Wavelength.
    dl_path : str, pathlib.Path
        Directory to search for files and download files to. If `None` defaults
        to the sunpy download directory.

    Returns
    -------
    sunpy.map.EUVIMap

    Notes
    -----
    This function tries to download a calibrated and archived file. If that
    fails it tries to download a beacon file.
    """
    if dl_path is None:
        dl_path = sunpy_config.get('downloads', 'download_dir')
    dl_path = Path(dl_path)

    map_path = _stereo_map_path(dtime, dl_path, wlen)
    if map_path.exists():
        return Map(map_path)

    try:
        dl_path = helpers.start_of_day_map(
            dtime, a.Instrument('EUVI'), a.Wavelength(195 * u.Angstrom))
    except Exception:
        map_path = _stereo_map_path(dtime, dl_path, wlen, beacon='_beacon')
        if map_path.exists():
            return Map(map_path)
        dl_path = stereo_beacon_start_of_day(dtime, dl_path, wlen)

    dl_path.replace(map_path)
    return Map(map_path)

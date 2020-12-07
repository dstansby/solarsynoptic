from datetime import datetime, timedelta

from sunpy.net import Fido, attrs as a
from sunpy.map import Map
from sunpy.time import parse_time

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
    if dtime > datetime.now():
        raise RuntimeError("Can't fetch map in the future!")
    # Get datetime at start of current day
    dtime = datetime(dtime.year, dtime.month, dtime.day)

    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(wlen))
    result = Fido.search(*query)
    try:
        download_path = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No maps available for whole day on {dtime}')

    return Map(download_path)

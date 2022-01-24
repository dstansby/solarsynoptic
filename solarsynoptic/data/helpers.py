from datetime import datetime, timedelta

from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time


def start_of_day(dtime):
    if isinstance(dtime, Time):
        dtime = dtime.datetime
    # Get datetime at start of current day
    return datetime(dtime.year, dtime.month, dtime.day)


def start_of_day_map(dtime, *query):
    """
    Download the first map available on a given date, with the given Fido query.

    Parameters
    ----------
    dtime :
        Datetime.
    query :
        Query parameters.
    """
    dtime = parse_time(dtime).to_datetime()
    if dtime > datetime.now():
        raise RuntimeError("Can't fetch map in the future!")
    # Get datetime at start of current day
    dtime = start_of_day(dtime)

    query = (a.Time(dtime, dtime + timedelta(days=1), dtime), *query)
    result = Fido.search(*query)
    try:
        download_path = Fido.fetch(result[0, 0])[0]
    except IndexError:
        raise RuntimeError(f'No maps available for whole day on {dtime}')

    return download_path

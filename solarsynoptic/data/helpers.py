from datetime import datetime, timedelta

from sunpy.net import Fido, attrs as a
from sunpy.map import Map
from sunpy.time import parse_time


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
    dtime = datetime(dtime.year, dtime.month, dtime.day)

    query = (a.Time(dtime, dtime + timedelta(days=1), dtime), *query)
    result = Fido.search(*query)
    try:
        download_path = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No maps available for whole day on {dtime}')

    return Map(download_path)

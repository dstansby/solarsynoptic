from datetime import datetime


def start_of_day(dtime):
    """
    Return a datetime at time 00:00:00 on the given date.
    """
    return datetime(dtime.year, dtime.month, dtime.day)

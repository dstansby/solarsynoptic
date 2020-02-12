from datetime import datetime


def start_of_day(dtime):
    return datetime(dtime.year, dtime.month, dtime.day)

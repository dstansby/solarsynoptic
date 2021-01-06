from sunpy.net import attrs as a

from . import helpers

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
    query = a.Instrument('AIA'), a.Wavelength(wlen)
    return helpers.start_of_day_map(dtime, *query)

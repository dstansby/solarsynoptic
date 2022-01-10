import astropy.units as u
import numpy as np
import sunpy.coordinates.sun
import sunpy.map

__all__ = ['weights_longitude']


def weights_longitude(width):
    """
    Weights to use when adding synoptic maps.

    Parameters
    ----------
    width : astropy.units.Quantity
        Width of the weight function.
    """

    def weights(smap):
        coord = sunpy.map.all_coordinates_from_map(smap)
        longs = coord.lon.to(u.deg).value
        l0 = smap.observer_coordinate.transform_to(smap.coordinate_frame).lon.to_value(u.deg)

        dcenterlong = (longs - l0 + 180) % 360 - 180
        weights = np.exp(-(dcenterlong / width.to_value(u.deg))**2)
        weights[weights < 0] = 0
        return weights / np.nanmax(weights)

    return weights

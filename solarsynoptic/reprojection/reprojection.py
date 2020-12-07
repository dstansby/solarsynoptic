"""
Tools for reprojecting maps.
"""
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from sunpy.map import Map, make_fitswcs_header
from reproject import reproject_interp

__all__ = ['reproject_carrington']


def reproject_carrington(smap, shape_out, latitude_projection='CAR'):
    """
    Reproject *smap* into a Carrington frame of reference.

    Notes
    -----
    The input map is reprojected into a map covering the full surface of the
    Sun, regardless of the extent of the input map.

    Parameters
    ----------
    smap : sunpy.map.GenericMap
    shape_out : [int, int]
        Number of output pixels in (latitude, longitude).

    Returns
    -------
    carrington_map : sunpy.map.Genericmap
        Reprojected map.
    """
    header_out = carrington_header(smap.date, shape_out,
                                   projection_code=latitude_projection)
    wcs_out = WCS(header_out)
    wcs_out.heliographic_observer = smap.observer_coordinate
    # Do the reprojection
    array, footprint = reproject_interp(smap, wcs_out, shape_out=shape_out)
    # Create a sunpy map, and copy over plot settings
    map_out = Map((array, header_out))
    map_out.plot_settings = smap.plot_settings
    return map_out


def carrington_header(dtime, shape_out, projection_code='CAR'):
    """
    Construct a FITS header for a Carrington coordinate frame.

    Parameters
    ----------
    dtime : astropy.time.Time
    shape_out : [int, int]
        Map shape in (latitude, longitude).
    projection_code : str, optional
        Projection to use for the latitude axis.

    Returns
    -------
    sunpy.map.MetaDict
    """
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime,
                         observer='earth')
    header = make_fitswcs_header(
        shape_out, frame_out,
        scale=[180 / shape_out[0],
               360 / shape_out[1]] * u.deg / u.pix,
        projection_code=projection_code)
    return header

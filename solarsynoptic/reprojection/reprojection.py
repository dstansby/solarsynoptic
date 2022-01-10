"""
Tools for reprojecting maps.
"""
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import numpy as np
import sunpy.coordinates
from sunpy.map import Map, make_fitswcs_header
from reproject import reproject_interp

__all__ = ['reproject_carrington']


def reproject_carrington(smap, shape_out, latitude_projection='CAR',
                         cache=True):
    """
    Reproject *smap* into a Carrington frame of reference.

    This is done using the `reproject.reproject_interp` function.

    Notes
    -----
    The input map is reprojected into a map covering the full surface of the
    Sun, regardless of the extent of the input map.

    *All* of the metadata in ``smap`` is copied to the output map, apart from
    the new WCS information. The ``BLANK`` keyword is also removed if present,
    as reprojection always results in floating point data.

    Parameters
    ----------
    smap : sunpy.map.GenericMap
    shape_out : [int, int]
        Number of output pixels in (latitude, longitude).
    cache : bool

    Returns
    -------
    carrington_map : sunpy.map.Genericmap
        Reprojected map.
    """
    if cache:
        from solarsynoptic.reprojection.database import DATABASE
        if smap in DATABASE:
            return DATABASE[smap]

    header_out = carrington_header(smap.date, shape_out,
                                   smap.observer_coordinate,
                                   projection_code=latitude_projection)
    wcs_out = WCS(header_out)
    # Do the reprojection
    array, footprint = reproject_interp(smap, wcs_out, shape_out=shape_out)
    # Copy some metadata over
    header_out['wavelnth'] = smap.meta.get('wavelnth', '')
    header_out['waveunit'] = smap.meta.get('waveunit', '')
    header_out['detector'] = smap.meta.get('detector', '')
    header_out['obsrvtry'] = smap.meta.get('obsrvtry', '')
    header_out['telescop'] = smap.meta.get('telescop', '')
    map_out = Map((array, header_out))
    map_out.plot_settings = smap.plot_settings

    if cache:
        from solarsynoptic.reprojection.database import save_and_cache
        save_and_cache(smap, map_out)
    return map_out


def carrington_header(dtime, shape_out, observer, projection_code='CAR'):
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
                         observer=observer)
    header = make_fitswcs_header(
        shape_out, frame_out,
        scale=[180 / shape_out[0],
               360 / shape_out[1]] * u.deg / u.pix,
        projection_code=projection_code)
    # Need to manually correct the CDELT2 value if using CEA
    if projection_code == 'CEA':
        # Since, this map uses the cylindrical equal-area (CEA) projection,
        # the spacing should be modified to 180/pi times the sin(latitude)
        # spacing
        # Reference: Section 5.5, Thompson 2006
        header['cdelt2'] = 180 / np.pi * 2 / shape_out[0]
    return header

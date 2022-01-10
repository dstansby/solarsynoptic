import astropy.units as u
import pytest
from astropy.coordinates import Longitude
from astropy.time import Time
from sunpy.map import GenericMap

from solarsynoptic.data import aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington


@pytest.mark.parametrize('latitude_projection', ['CEA', 'CAR'])
@pytest.mark.parametrize('cache', [True, False])
def test_reprojection(cache, latitude_projection):
    smap = aia_start_of_day_map(Time.now(), 193 * u.Angstrom)

    shape_out = (180, 360)
    car_map = reproject_carrington(smap, shape_out, cache=cache,
                                   latitude_projection=latitude_projection)

    assert isinstance(car_map, GenericMap)
    assert car_map.data.shape == shape_out
    ll_coord = car_map.pixel_to_world(*(-0.5, -0.5) * u.pix)
    assert u.allclose(ll_coord.lon, Longitude(-180 * u.deg))
    assert u.allclose(ll_coord.lat, -90 * u.deg)

    assert car_map.meta['CTYPE1'][5:] == latitude_projection
    assert car_map.meta['CTYPE2'][5:] == latitude_projection

    assert car_map.plot_settings == smap.plot_settings

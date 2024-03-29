import astropy.units as u
import pytest
from astropy.time import Time
from sunpy.map.sources import AIAMap

from solarsynoptic.data.aia import aia_start_of_day_map


@pytest.mark.parametrize('t', [Time('2020-01-01'), Time.now()])
def test_aia_start_of_day_map(t):
    if t < Time('2021-01-01'):
        pytest.xfail('Older AIA data has server issues')
    smap = aia_start_of_day_map(t, 193 * u.Angstrom)
    assert isinstance(smap, AIAMap)

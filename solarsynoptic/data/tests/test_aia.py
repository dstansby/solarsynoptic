import pytest

import astropy.units as u
from astropy.time import Time
from sunpy.map.sources import AIAMap

from solarsynoptic.data.aia import aia_start_of_day_map


@pytest.mark.parametrize('t', [Time('2020-01-01'), Time.now()])
def test_aia_start_of_day_map(t):
    smap = aia_start_of_day_map(t, 193 * u.Angstrom)
    assert isinstance(smap, AIAMap)

import astropy.units as u
import pytest
from astropy.time import Time
from sunpy.map.sources import EUVIMap

from solarsynoptic.data import stereo_start_of_day_map


@pytest.mark.parametrize('t', [Time('2020-01-01'), Time.now() - 2*u.day])
def test_stereo_start_of_day_map(t):
    smap = stereo_start_of_day_map(t)
    assert isinstance(smap, EUVIMap)

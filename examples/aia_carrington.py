"""
AIA Carrington map
==================
"""

import astropy.units as u

from solarsynoptic.data import aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington

shape_out = [180, 360]
map_in = aia_start_of_day_map('2020-01-01', 193 * u.Angstrom)
map_out = reproject_carrington(map_in, shape_out)
map_out.peek()

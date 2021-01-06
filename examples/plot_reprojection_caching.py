"""
AIA Carrington map
==================
"""
###############################################################################
# Import required functions
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map

from solarsynoptic.data import aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington
from solarsynoptic.reprojection.database import DATABASE

###############################################################################
# Define the output of the shape map. This is number of pixels in latitude and
# longitude.
shape_out = [180, 360]
###############################################################################
# Fetch a single AIA map. This will download the first map observed on the
# given day.
# map_in = aia_start_of_day_map('2020-01-01', 193 * u.Angstrom)
map_in = sunpy.map.Map('/Users/dstansby/sunpy/data/aia_lev1_193a_2020_01_01t00_00_04_84z_image_lev1.fits')
map_out = reproject_carrington(map_in, [180, 360], cache=True)
print(map_out)

"""
AIA Carrington map
==================
"""
###############################################################################
# Import required functions
import astropy.units as u
import matplotlib.pyplot as plt

from solarsynoptic.data import aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington

###############################################################################
# Define the output of the shape map. This is number of pixels in latitude and
# longitude.
shape_out = [180, 360]
###############################################################################
# Fetch a single AIA map. This will download the first map observed on the
# given day.
map_in = aia_start_of_day_map('2020-01-01', 193 * u.Angstrom)
###############################################################################
# Reproject the map into a Carrington frame of reference, and show the result.
map_out = reproject_carrington(map_in, shape_out)
fig = plt.figure()
map_out.plot()

###############################################################################
# By default `reproject_carrington` uses a Carée projection, where the latitude
# is evenly sampled. We can change this to an equal area (CEA) projection,
# where sin(latitude) is evenly sampled instead.
map_out = reproject_carrington(map_in, shape_out, latitude_projection='CEA')
fig = plt.figure()
map_out.plot()

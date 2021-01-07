"""
AIA Carrington map
==================
"""
###############################################################################
# Import required functions
from datetime import datetime, timedelta

import astropy.units as u
import matplotlib.pyplot as plt

from solarsynoptic.data import aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington
from solarsynoptic.coadd import coadd

###############################################################################
# Define the output of the shape map. This is number of pixels in latitude and
# longitude.
shape_out = [180, 360]
###############################################################################
# Fetch a single AIA map. This will download the first map observed on the
# given day.
maps_in = []
for i in range(24):
    d = datetime.now() - timedelta(days=i)
    maps_in.append(aia_start_of_day_map(d, 193 * u.Angstrom))

###############################################################################
# Reproject the maps into a Carrington frame of reference
maps_in = [reproject_carrington(map_in, shape_out, latitude_projection='CEA') for map_in in maps_in]

###############################################################################
# Add the maps together
map_out = coadd(maps_in)
fig = plt.figure()
map_out.plot()

plt.show()

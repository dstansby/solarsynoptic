"""
AIA Carrington map
==================
"""
###############################################################################
# Import required functions
from datetime import datetime, timedelta
import numpy as np

import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map

from solarsynoptic.data import stereo_start_of_day_map, aia_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington
from solarsynoptic.coadd import coadd, long_weights

###############################################################################
# Define the output of the shape map. This is number of pixels in latitude and
# longitude.
shape_out = [180, 360]

###############################################################################
# Download STEREO maps
maps_in_stereo = []
for i in range(10):
    d = datetime.now() - timedelta(days=i)
    smap = stereo_start_of_day_map(d)
    # Normalise to AIA
    data = (smap.data - 720.0) * 110 / 170
    data[data < 0] = np.nan
    smap = sunpy.map.Map((data, smap.meta))
    maps_in_stereo.append(smap)

###############################################################################
# Download AIA maps
maps_in_aia = []
for i in range(17):
    d = datetime.now() - timedelta(days=i)
    maps_in_aia.append(aia_start_of_day_map(d, 193 * u.Angstrom))

norm = maps_in_aia[0].plot_settings['norm']
norm.vmin = 40
norm.vmax = 5000

###############################################################################
# Reproject the maps into a Carrington frame of reference
maps_in_stereo = [reproject_carrington(map_in, shape_out) for map_in in maps_in_stereo]
maps_in_aia = [reproject_carrington(map_in, shape_out) for map_in in maps_in_aia]


def weight_function(smap):
    weights = long_weights(30 * u.deg)(smap)
    factor = (datetime.now() - smap.date.to_datetime()) / timedelta(days=1)
    return weights / factor

###############################################################################
# Add the maps together
map_out = coadd(maps_in_stereo, weight_function=weight_function)
fig = plt.figure()
map_out.plot(cmap='sdoaia193', norm=norm)
fig.savefig('figs/stereo.png')

map_out = coadd(maps_in_aia, weight_function=weight_function)
fig = plt.figure()
map_out.plot(cmap='sdoaia193', norm=norm)
fig.savefig('figs/aia.png')

map_out = coadd(maps_in_aia + maps_in_stereo, weight_function=weight_function)
fig = plt.figure()
map_out.plot(cmap='sdoaia193', norm=norm)
fig.savefig('figs/combined.png')

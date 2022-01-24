"""
AIA Carrington map
==================
"""
###############################################################################
# Import required functions
from datetime import datetime, timedelta

import astropy.units as u
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map

from solarsynoptic.combine import coadd, weights_longitude
from solarsynoptic.data import aia_start_of_day_map, stereo_start_of_day_map
from solarsynoptic.reprojection import reproject_carrington

###############################################################################
# Define the output of the shape map. This is number of pixels in latitude and
# longitude.
shape_out = [720, 1440]

###############################################################################
# Download STEREO maps
maps_in_stereo = []
for i in range(2):
    d = datetime.now() - timedelta(days=i) - timedelta(days=1)
    smap = stereo_start_of_day_map(d)
    # Normalise to AIA
    data = (smap.data - 720.0) * 110 / 170
    data[data < 0] = np.nan
    smap = sunpy.map.Map((data, smap.meta))
    maps_in_stereo.append(smap)

###############################################################################
# Download AIA maps
maps_in_aia = []
for i in range(2):
    d = datetime.now() - timedelta(days=i) - timedelta(days=1)
    maps_in_aia.append(aia_start_of_day_map(d, 193 * u.Angstrom))

norm = maps_in_aia[0].plot_settings['norm']
norm.vmin = 40
norm.vmax = 5000

###############################################################################
# Reproject the maps into a Carrington frame of reference
maps_in_stereo = [reproject_carrington(map_in, shape_out) for map_in in maps_in_stereo]
maps_in_aia = [reproject_carrington(map_in, shape_out) for map_in in maps_in_aia]


def weight_function(smap):
    weights = weights_longitude(30 * u.deg)(smap)
    factor = (datetime.now() - smap.date.to_datetime()) / timedelta(days=1)
    return weights / factor


###############################################################################
# Add the maps together
dtime = datetime.now().strftime('%Y-%m-%d')

# Just the STEREO maps
map_out = coadd(maps_in_stereo, weight_function=weight_function)
fig = plt.figure()
map_out.plot(cmap='sdoaia193', norm=norm)
plt.gca().set_title(f'EUVI, updated {dtime}')

# Just tthe AIA maps
map_out = coadd(maps_in_aia, weight_function=weight_function)
fig = plt.figure()
map_out.plot(cmap='sdoaia193', norm=norm)
map_out.save(f'maps/aia_synoptic_{dtime}.fits')
plt.gca().set_title(f'AIA, updated {dtime}')

# STEREO and AIA maps
map_out = coadd(maps_in_aia + maps_in_stereo, weight_function=weight_function)
map_out.save(f'maps/aia_euvi_synoptic_{dtime}.fits', overwrite=True)

fig = plt.figure()
map_out = sunpy.map.Map((norm(map_out.data), map_out.meta))
map_out.plot(cmap='sdoaia193', norm=mcolor.Normalize(vmin=0, vmax=1))
plt.gca().set_title(f'AIA + EUVI, udpated {dtime}')

plt.show()

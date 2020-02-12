from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.map import Map
# Register sunpy colourmaps
import sunpy.visualization.colormaps

from aia_helpers import create_synoptic_map


figdir = Path('synoptic_maps')


def save_map(m, crot):
    m.save(figdir / 'fits' / f'aia193_synoptic_{crot}.fits', overwrite=True)


def save_figure(m, crot):
    fig = plt.figure()
    m.plot(cmap='sdoaia193')
    # Plot formatting
    ax = plt.gca()
    # Add start and end times to title
    tstart = carrington_rotation_time(crot - 1).to_datetime().strftime(
        '%Y-%m-%d')
    tend = carrington_rotation_time(crot).to_datetime().strftime(
        '%Y-%m-%d')
    ax.set_title(f'Carrington rotation {crot} ({tstart} - {tend})')
    # This has to be slightly larger than -0.5 for the y-axis labels to show
    ax.set_xlim(left=-0.49999)
    ax.grid(False)
    # Save plot
    plt.savefig(figdir / 'plots' / f'aia193_synoptic_{crot}.pdf',
                bbox_inches='tight')
    plt.close('all')


if __name__ == '__main__':
    crot = 2181

    while True:
        t = carrington_rotation_time(crot).to_datetime()
        map = create_synoptic_map(t)
        # Norm the data
        data = map.data
        data = map.plot_settings['norm'](data)
        # Save the map
        m = Map((data, map.meta))
        save_map(m, crot)
        save_figure(m, crot)

        crot -= 1

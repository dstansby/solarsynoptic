import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.map import Map
# Register sunpy colourmaps
import sunpy.visualization.colormaps

from aia_helpers import create_synoptic_map

# The directory to save fits files and figures to
output_dir = Path('synoptic_maps')


def save_map(m, crot):
    m.save(output_dir / 'fits' / f'aia193_synoptic_{crot}.fits', overwrite=True)


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
    plt.savefig(output_dir / 'plots' / f'aia193_synoptic_{crot}.pdf',
                bbox_inches='tight')
    plt.close('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create AIA Carrington maps.')
    parser.add_argument('--crot', metavar='crot', type=int, nargs=1,
                        help='Carrington rotation to start at '
                             '(before working backwards).', required=True)
    parser.add_argument('--wlen', metavar='wavelegnth', type=int, nargs=1,
                        help='Wavelength to use. Must be 193,', required=True)
    args = parser.parse_args()
    crot = args.crot[0]
    wlen = args.wlen[0]

    wlens = [193, 211]
    assert wlen in wlens, f'Wavelength must be in {wlens}'

    while True:
        t = carrington_rotation_time(crot).to_datetime()
        map = create_synoptic_map(t, wlen)
        # Norm the data
        data = map.data
        data = map.plot_settings['norm'](data)
        # Save the map
        m = Map((data, map.meta))
        save_map(m, crot)
        save_figure(m, crot)

        crot -= 1

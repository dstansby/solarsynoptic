import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.map import Map
# Register sunpy colourmaps
import sunpy.visualization.colormaps

from aia_helpers import create_synoptic_map

# The directory to save fits files and figures to
output_dir = Path('/Volumes/Work/Data/synoptic_aia/synoptic_fits')
(output_dir / 'fits').mkdir(exist_ok=True)
(output_dir / 'plots').mkdir(exist_ok=True)


def save_map(m, crot, wlen):
    m.save(output_dir / 'fits' / f'aia{wlen}_synoptic_{crot}.fits',
           overwrite=True)


def save_figure(m, crot, wlen):
    fig = plt.figure()
    m.plot(cmap=f'sdoaia{wlen}', vmax=2000)
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
    plt.savefig(output_dir / 'png' / f'aia{wlen}_synoptic_{crot}.png', dpi=200)
    plt.close('all')


if __name__ == '__main__':
    wlens = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]
    parser = argparse.ArgumentParser(description='Create AIA Carrington maps.')
    parser.add_argument('--crot', metavar='crot', type=int, nargs=1,
                        help='Carrington rotation to start at '
                             '(before working backwards).', required=True)
    parser.add_argument('--wlen', metavar='wavelegnth', type=int, nargs=1,
                        help=f'Wavelength to use. Must be in {wlens},',
                        required=True)
    args = parser.parse_args()
    crot = args.crot[0]
    wlen = args.wlen[0]

    assert wlen in wlens, f'Wavelength must be in {wlens}'

    while True:
        t = carrington_rotation_time(crot).to_datetime()
        m = create_synoptic_map(t, wlen)
        print(f'Created AIA {wlen} synoptic map for rotation {crot} 🎉')
        # Save the map
        save_map(m, crot, wlen)
        save_figure(m, crot, wlen)

        crot += 1

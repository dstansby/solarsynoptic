Getting started
---------------
When producing synoptic maps of the Sun, there are typically 3 steps involved:

- Source the desired images, in their original Helioprojective coordinate system.
- Reproject each image into a map projection that spans the full solar surface.
- Add the images together to cover a large fractin of the solar surface

solarsynoptic breaks these tasks down into three sub-modules:

- `solarsynoptic.data`
- `solarsynoptic.reprojection`
- `solarsynoptic.coadd`

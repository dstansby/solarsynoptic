=============
solarsynoptic
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   auto_examples/index


Getting started
===============
When producing synoptic maps of the Sun, there are typically 3 steps involved:

- Source the desired images, in their original Helioprojective coordinate system.
- Reproject each image into a map projection that spans the full solar surface.
- Add the images together to cover a large fraction of the solar surface

solarsynoptic breaks these tasks down into three sub-modules:

- `solarsynoptic.data`
- `solarsynoptic.reprojection`
- `solarsynoptic.combine`


Fetching data
-------------
`solarsynoptic.data` contains a couple of helper functions for downloading
AIA and STEREO maps from the beginning of a given day. For recent dates they
automatically handle getting a near-real-time image when a final archived
image isn't yet available.

.. automodapi:: solarsynoptic.data
  :no-heading:

Reprojecting data
-----------------
`solarsynoptic.reprojection` contains a single function to reproject maps into
a Carrington coordinate frame.

.. automodapi:: solarsynoptic.reprojection
  :no-heading:

Caching reprojected maps
~~~~~~~~~~~~~~~~~~~~~~~~
Under the hood ``solarsynoptic`` can cache the reprojected map to a file on
disk, which avoids having to run the reprojection routine each time the
reprojected map is requested. This is enabled by default, and can be disabled
by passing ``cache=False`` to ``reproject_carrington``.

Combining multiple maps
-----------------------

.. automodapi:: solarsynoptic.combine
  :no-heading:

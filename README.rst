Solar Synoptic Maps
===================

Code to create synoptic maps from Solar Dynamics Observatory (SDO) data.

Getting started
---------------
Before starting, you should set:
  - The directory in which maps are downloaded to and read from at the top of ``aia_helpers.py``
  - The directory to save fits files and figures to at the top of ``aia_synoptic_maps.py``

Running
-------
Example run:

```
python aia_synoptic_maps.py --crot 2146 --wlen 193
```

Allowed wavelengths are ``[94, 131, 171, 193, 211, 304, 335, 1600, 1700]``.

Because reprojecting the maps is relatively expensive, Carrington projections
of the maps are automatically saved under ``aia_193_synoptic_{date}.fits``. This
makes it fast to modify the method by which multiple maps are added together
to form the final synoptic map.

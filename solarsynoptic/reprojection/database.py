"""
Tools for caching reprojected images.
"""
import json
import hashlib
from pathlib import Path

import astropy.units as u
import sunpy.util.config
import pandas as pd


def dict_hash(dictionary):
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # We need to sort arguments so {'a': 1, 'b': 2} is
    # the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


def item_hash(smap):
    """
    Create a hash based on the stored items in a map MetaDict.

    This relies on all the items themselves being hashable. For this reason
    the 'keycomments' item, which is a dict, is excluded from the hash.

    If creating the hash fails, returns `None`.
    """
    meta = smap.meta
    meta_copy = meta.copy()
    meta_copy.pop('keycomments', None)
    return dict_hash(meta_copy)


class _DataBase:
    def __init__(self):
        # Get the directory
        sunpy_dl_dir = sunpy.util.config.get_and_create_download_dir()
        solarsynop_dir = Path(sunpy_dl_dir) / '..' / 'solarsynoptic'
        # Make some directories
        solarsynop_dir.mkdir(exist_ok=True)
        (solarsynop_dir / 'data').mkdir(exist_ok=True)
        self.directory = solarsynop_dir.resolve()

        if not self.db_fpath.exists():
            self.create_empty_db()

        self.db = self.load_db()

    @property
    def db_fpath(self):
        return self.directory / 'solarsynoptic_db.csv'

    def create_empty_db(self):
        db = pd.DataFrame(columns=['Original file hash',
                                   'Reprojected file hash',
                                   'Reprojected file path'],
                          dtype=str)
        db = db.set_index('Original file hash')
        db.to_csv(self.db_fpath)

    def load_db(self):
        """
        Load the solarsynoptic database file.

        Returns
        -------
        pandas.DataFrame
        """
        db = pd.read_csv(self.db_fpath, dtype=str)
        db = db.set_index('Original file hash')
        return db

    def save_db(self):
        """
        Save the database to file.
        """
        self.db.to_csv(self.db_fpath)

    def insert(self, original_map, projection,
               reprojected_map, reprojected_filepath):
        # Add entry to database
        original_meta_hash = self.key(original_map, projection)
        reprojected_meta_hash = item_hash(reprojected_map)
        self.db.loc[original_meta_hash] = [reprojected_meta_hash,
                                           str(reprojected_filepath)]
        self.save_db()

    def __contains__(self, item):
        return self.key(*item) in self.db.index

    def __getitem__(self, item):
        fpath = self.db.loc[self.key(*item), 'Reprojected file path']
        return sunpy.map.Map(fpath)

    @staticmethod
    def key(smap, projection):
        # Return database key for a given map and output projection
        return item_hash(smap) + projection


DATABASE = _DataBase()


def save_and_cache(original_map, projection, reprojected_map):
    """
    Save ``reprojected_map`` to a .FITS file, and insert an entry into the
    solarsynoptic database.

    Parameters
    ----------
    original_map : sunpy.map.GenericMap
        The original (pre-reprojection) map.
    reprojected_map : sunpy.map.GenericMap
        The reprojected map.
    """
    # Save file
    ctype = reprojected_map.meta['CTYPE2'].upper()
    reprojected_fname = (construct_fname(original_map) +
                         f'_reprojected_{ctype}.fits')
    reprojected_file_path = DATABASE.directory / 'data' / reprojected_fname
    reprojected_map.save(reprojected_file_path, overwrite=True)

    # Add to database
    DATABASE.insert(original_map, projection,
                    reprojected_map, reprojected_file_path)


def construct_fname(smap):
    """
    Construct a human readable and (hopefully!) unique filename from a map.

    This takes the form {observatory}_{detector}_{wavelength}_{timestamp}.

    Parameters
    ----------
    smap : sunpy.map.GenericMap

    Returns
    -------
    fname : str
        Filename. Does not include any filename extension.
    """
    datestr = smap.date.strftime('%Y%m%d-%H%M%S-%f')
    wlen = int(smap.measurement.to_value(u.Angstrom))
    return (f'{smap.observatory}_{smap.detector}_'
            f'{wlen}_{datestr}')

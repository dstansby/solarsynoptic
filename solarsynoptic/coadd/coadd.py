import numpy as np
import sunpy.map


def coadd(smaps, weight_function=None):
    r"""
    Add a set of full Sun maps together.

    This is helpful for taking several observations, and adding them together
    to increase coverage of the Sun.

    Parameters
    ----------
    smaps : list[sunpy.map.GenericMap]
    weight_function : callable
        A function that takes a single `~sunpy.map.GenericMap` as input, and
        returns an array of weights that is the same shape as the input map.
    """
    smaps = sunpy.map.MapSequence(*smaps)
    if not smaps.all_maps_same_shape():
        raise ValueError('Not all the input maps are the same shape.')

    out_shape = smaps[0].data.shape
    out_data = np.zeros(out_shape)
    out_weights = np.zeros(out_shape)
    # Mask to keep track of NaNs. Is set to 0 if all summed pixels are NaN,
    # 1 if at least one summed pixel has a none NaN value.
    non_nan_mask = np.zeros(out_shape)
    for smap in smaps:
        if weight_function is None:
            weights = np.isfinite(smap.data)
        else:
            weights = weight_function(smap)
            # Validate the weights
            if weights.shape != out_shap:
                raise RuntimeError(
                    f'weight_function returned data shape {weights.shape}, '
                    f'but expected shape is {out_shape}')
            if np.any(~np.isfinite(weights)):
                raise RuntimeError(
                    'weight_funciton returned at least one non-finite value')

        non_nan_mask = np.logical_or(non_nan_mask, ~np.isnan(smap.data))
        out_data += np.nan_to_num(smap.data) * weights
        out_weights += weights

    # Normalise to the total weights
    # Ignore divide by zero
    with np.errstate(divide='ignore'):
        out_data /= out_weights
    out_data[~non_nan_mask] = np.nan

    out_map = sunpy.map.Map((out_data, smaps[-1].meta))
    return out_map

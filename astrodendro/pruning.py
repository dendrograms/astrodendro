# Licensed under an MIT open source license - see LICENSE
"""
The pruning module provides several functions to perform
common pruning via the ``is_independent`` keyword in the Dendrogram :meth:`~astrodendro.dendrogram.Dendrogram.compute` method.

Examples::

    #prune unless leaf peak value >= 5
    Dendrogram.compute(data, is_independent=min_peak(5))

    #prune unless leaf contains 10 pixels
    Dendrogram.compute(data, is_independent=min_npix(10))

    #apply both criteria
    is_independent = all_true((min_peak(5), min_npix(10)))
    Dendrogram.compute(data, is_independent=is_independent)
"""

import numpy as np

def all_true(funcs):
    """Combine several ``is_independent`` functions into one

    Parameters
    ----------
    funcs : list-like
        A list of ``is_independent`` functions

    Returns
    -------
    combined_func : function
        A new function which returns true of all the input functions are true
    """
    def result(*args, **kwargs):
        return all(f(*args, **kwargs) for f in funcs)
    return result


def min_delta(delta):
    """
    Minimum delta criteria

    Parameters
    ----------
    delta : float
        The minimum height of a leaf above its merger level

    """
    def result(structure, index=None, value=None):
        if value is None:
            return (structure.vmax - structure.vmin) >= delta
        return (structure.vmax - value) >= delta
    return result


def min_sum(sum):
    """
    Minimum sum criteria

    Parameters
    ----------
    sum : float
        The minimum sum of the pixel values in a leaf
    """
    def result(structure, index=None, value=None):
        return np.nansum(structure.values()) >= sum
    return result


def min_peak(peak):
    """
    Minimum peak criteria

    Parameters
    ----------
    peak : float
        The minimum peak pixel value in a leaf
    """
    def result(structure, index=None, value=None):
        return structure.vmax >= peak
    return result


def min_npix(npix):
    """
    Minimum npix criteria

    Parameters
    ----------
    npix : int
        The minimum number of pixels in a leaf
    """
    def result(structure, index=None, value=None):
        return len(structure.values()) >= npix
    return result


def contains_seeds(seeds):
    """
    Critieria that leaves contain at least one of a list of seed positions

    Parameters
    ----------
    seeds : tuple of array-like
        seed locations. The ith array in the tuple lists the ith coordinate
        for each seed. This is the format returned, e.g., by np.where
    """
    shp = [np.asarray(s).max() + 2 for s in seeds]
    rav = np.ravel_multi_index(seeds, shp)

    def result(structure, index=None, value=None):
        sid = structure.indices()
        if len(sid) != len(seeds):
            raise TypeError("Dimensions of seeds and data do not agree")
        rav2 = np.ravel_multi_index(sid, shp, mode='clip')
        return np.intersect1d(rav, rav2).size > 0

    return result

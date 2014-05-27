
import numpy as np
from . import pruning
from . import six
from .dendrogram import _sorted_by_idx
from copy import copy

def post_pruning(d, min_value=-np.inf, min_delta=0, min_npix=0,
                is_independent=None):
    '''
    Prune a dendrogram after it has been computed.

    Parameters
    ----------
    d : Dendrogram
        Computed dendrogram object.
    min_value : float, optional
        The minimum data value to go down to when computing the
        dendrogram. Values below this threshold will be ignored.
    min_delta : float, optional
        The minimum height a leaf has to have in order to be considered an
        independent entity.
    min_npix : int, optional
        The minimum number of pixels/values needed for a leaf to be considered
        an independent entity.
    is_independent : function or list of functions, optional
        A custom function that can be specified that will determine if a
        leaf can be treated as an independent entity. The signature of the
        function should be ``func(structure, index=None, value=None)``
        where ``structure`` is the structure under consideration, and
        ``index`` and ``value`` are optionally the pixel that is causing
        the structure to be considered for merging into/attaching to the
        tree.

        If multiple functions are provided as a list, they
        are all applied when testing for independence.
    '''

    d_new = copy(d)

    tests = [pruning.min_delta(min_delta),
             pruning.min_npix(min_npix)]
    if is_independent is not None:
        if hasattr(is_independent, '__iter__'):
            tests.extend(is_independent)
        else:
            tests.append(is_independent)
    is_independent = pruning.all_true(tests)

    keep_structures = {}
    for struct in d.all_structures:
        if is_independent(struct):
            keep_structures[struct.idx] = struct


            # if struct.is_branch:
            #     children = struct.children
            #     if not children:
            #         struct.is_leaf = True

    # Create trunk from objects with no ancestors
    d_new.trunk = _sorted_by_idx([structure for structure in six.itervalues(keep_structures) if structure.parent is None])
    print d_new.trunk
    # Remove orphan leaves that aren't large enough
    leaves_in_trunk = [structure for structure in d_new.trunk if structure.is_leaf]
    for leaf in leaves_in_trunk:
        if not is_independent(leaf):
            # This leaf is an orphan, so remove all references to it:
            structures.pop(leaf.idx)
            d_new.trunk.remove(leaf)
            leaf._fill_footprint(d_new.index_map, -1)

    # To make the structure.level property fast, we ensure all the structures in the
    # trunk have their level cached as "0"
    for structure in d_new.trunk:
        structure._level = 0  # See the definition of level() in structure.py

    # Save a list of all structures accessible by ID
    d_new._structures_dict = {}

    # Re-assign idx and update index map
    sorted_structures = sorted(d, key=lambda s: s.smallest_index)
    for idx, s in enumerate(sorted_structures):
        s.idx = idx
        s._fill_footprint(d_new.index_map, idx, recursive=False)
        d_new._structures_dict[idx] = s

    print d_new.__len__()

    return d_new

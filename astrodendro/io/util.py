import numpy as np

from .. import six
from astropy.utils.console import ProgressBar
from astropy import log
import time


def parse_newick(string):

    items = {}

    # Find maximum level
    current_level = 0
    max_level = 0
    log.debug('String loading...')
    for i, c in enumerate(string):
        if c == '(':
            current_level += 1
        if c == ')':
            current_level -= 1
        max_level = max(max_level, current_level)

    # Loop through levels and construct tree
    log.debug('Tree loading...')
    for level in range(max_level, 0, -1):

        pairs = []

        current_level = 0
        for i, c in enumerate(string):
            if c == '(':
                current_level += 1
                if current_level == level:
                    start = i
            if c == ')':
                if current_level == level:
                    pairs.append((start, i))
                current_level -= 1

        for pair in pairs[::-1]:

            # Extract start and end of branch definition
            start, end = pair

            # Find the ID of the branch
            colon = string.find(":", end)
            branch_id = string[end + 1:colon]
            if branch_id == '':
                branch_id = 'trunk'
            else:
                branch_id = int(branch_id)

            # Add branch definition to overall definition
            items[branch_id] = eval("{%s}" % string[start + 1:end])

            # Remove branch definition from string
            string = string[:start] + string[end + 1:]

    new_items = {}

    def collect(d):
        for item in d:
            if item in items:
                collect(items[item])
                d[item] = (items[item], d[item])
        return

    collect(items['trunk'])

    return items['trunk']


def parse_dendrogram(newick, data, index_map, readmethod=3):
    from ..dendrogram import Dendrogram
    from ..structure import Structure

    t0 = time.time()

    d = Dendrogram()
    d.ndim = len(data.shape)

    d._structures_dict = {}
    d.data = data
    d.index_map = index_map

    flux_by_structure = {}
    indices_by_structure = {}

    def _construct_tree(repr):
        structures = []
        for idx in repr:
            idx = int(idx)
            structure_indices = indices_by_structure[idx]
            f = flux_by_structure[idx]
            if type(repr[idx]) == tuple:
                sub_structures_repr = repr[idx][0]  # Parsed representation of sub structures
                sub_structures = _construct_tree(sub_structures_repr)
                for i in sub_structures:
                    d._structures_dict[i.idx] = i
                b = Structure(structure_indices, f, children=sub_structures, idx=idx, dendrogram=d)
                # Correct merge levels - complicated because of the
                # order in which we are building the tree.
                # What we do is look at the heights of this branch's
                # 1st child as stored in the newick representation, and then
                # work backwards to compute the merge level of this branch
                first_child_repr = six.next(six.itervalues(sub_structures_repr))
                if type(first_child_repr) == tuple:
                    height = first_child_repr[1]
                else:
                    height = first_child_repr
                d._structures_dict[idx] = b
                structures.append(b)
            else:
                l = Structure(structure_indices, f, idx=idx, dendrogram=d)
                structures.append(l)
                d._structures_dict[idx] = l
        return structures

    if readmethod==3:
        from scipy import ndimage
        """
        In [14]: %timeit np.unique(index_map)
        1 loops, best of 3: 1.81 s per loop

        In [15]: %timeit np.unique(index_map[index_map>-1])
        10 loops, best of 3: 105 ms per loop
        """
        idxs = np.unique(d.index_map[d.index_map > -1])

        # ndimage ignores 0 and -1, but we want index 0
        object_slices = ndimage.find_objects(d.index_map+1)
        index_cube = np.indices(d.index_map.shape)

        # Need to have same length, otherwise assumptions above are wrong
        assert len(idxs) == len(object_slices)

        for idx,sl in ProgressBar(zip(idxs, object_slices)):
            match = d.index_map[sl] == idx
            sl2 = (slice(None),) + sl
            match_inds = index_cube[sl2][:, match]
            coords = list(zip(*match_inds))
            data = d.data[sl][match].tolist()
            if idx in flux_by_structure:
                flux_by_structure[idx] += data
                indices_by_structure[idx] += coords
            else:
                flux_by_structure[idx] = data
                indices_by_structure[idx] = coords


    if readmethod==2:
        # Alternative implementation.  Turns out to be slower.
        indices = np.unique(d.index_map[d.index_map>-1])
        log.debug('[np way] Creating index maps for {0} indices...'.format(len(indices)))
        for idx in ProgressBar(indices):
            match = d.index_map == idx # This is probably why it's slower
            whmatch = np.nonzero(match)
            if idx in flux_by_structure:
                flux_by_structure[idx] += d.data[whmatch].tolist()
                indices_by_structure[idx] += zip(*whmatch)
            else:
                flux_by_structure[idx] = d.data[whmatch].tolist()
                indices_by_structure[idx] = zip(*whmatch)

    if readmethod == 1:
        # Do a fast iteration through d.data, adding the indices and data values
        # to the two dictionaries declared above:
        indices = np.array(np.where(d.index_map > -1)).transpose()

        log.debug('Creating index maps for {0} indices...'.format(len(indices)))
        for coord in ProgressBar(indices):
            coord = tuple(coord)
            idx = d.index_map[coord]
            if idx in flux_by_structure:
                flux_by_structure[idx].append(d.data[coord])
                indices_by_structure[idx].append(coord)
            else:
                flux_by_structure[idx] = [d.data[coord]]
                indices_by_structure[idx] = [coord]

    """
    In [7]: util.parse_dendrogram(newick, data, index_map, readmethod=1)
    |=========================================================================================================================| 330k/330k (100.00%)        28s
    Out[7]: <astrodendro.dendrogram.Dendrogram at 0x107678450>

    In [8]: util.parse_dendrogram(newick, data, index_map, readmethod=2)
    |=========================================================================================================================| 528 /528  (100.00%)      3m15s
    Out[8]: <astrodendro.dendrogram.Dendrogram at 0x1522e5f90>
    """

            

    log.debug('Parsing newick and constructing tree...')
    d.trunk = _construct_tree(parse_newick(newick))
    # To make the structure.level property fast, we ensure all the items in the
    # trunk have their level cached as "0"
    for structure in d.trunk:
        structure._level = 0  # See the @property level() definition in structure.py

    d._index()
    log.info("Read method {0} completed in {1} seconds.".format(readmethod,
                                                                time.time()-t0))
    return d

"""
Performance stats for a 528-element dendrogram:
INFO: Read method 1 completed in 51.4208781719 seconds. [astrodendro.io.util]
INFO: Read method 2 completed in 220.598875999 seconds. [astrodendro.io.util]
INFO: Read method 3 completed in 27.0453000069 seconds. [astrodendro.io.util]

For a 4000 element dendrogram:
INFO: Read method 3 completed in 30.4590411186 seconds. [astrodendro.io.util]
"""

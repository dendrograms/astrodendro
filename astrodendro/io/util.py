import numpy as np

from .. import six
from astropy.utils.console import ProgressBar
from astropy import log


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


def parse_dendrogram(newick, data, index_map):
    from ..dendrogram import Dendrogram
    from ..structure import Structure

    d = Dendrogram()
    d.ndim = len(data.shape)

    d._structures_dict = {}
    d.data = data
    d.index_map = index_map

    try:
        flux_by_structure, indices_by_structure = _fast_reader(d.index_map, data)
    except ImportError:
        flux_by_structure, indices_by_structure = _slow_reader(d.index_map, data)

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

    log.debug('Parsing newick and constructing tree...')
    d.trunk = _construct_tree(parse_newick(newick))
    # To make the structure.level property fast, we ensure all the items in the
    # trunk have their level cached as "0"
    for structure in d.trunk:
        structure._level = 0  # See the @property level() definition in structure.py

    d._index()
    return d

def _fast_reader(index_map, data):
    """
    Use scipy.ndimage.find_objects to quickly identify subsets of the data
    to increase speed of dendrogram loading
    """

    flux_by_structure, indices_by_structure = {},{}

    from scipy import ndimage
    idxs = np.unique(index_map[index_map > -1])

    # ndimage ignores 0 and -1, but we want index 0
    object_slices = ndimage.find_objects(index_map+1)
    index_cube = np.indices(index_map.shape)

    # Need to have same length, otherwise assumptions above are wrong
    assert len(idxs) == len(object_slices)
    log.debug('Creating index maps for {0} indices...'.format(len(idxs)))

    for idx,sl in ProgressBar(zip(idxs, object_slices)):
        match = index_map[sl] == idx
        sl2 = (slice(None),) + sl
        match_inds = index_cube[sl2][:, match]
        coords = list(zip(*match_inds))
        dd = data[sl][match].tolist()
        flux_by_structure[idx] = dd
        indices_by_structure[idx] = coords

    return flux_by_structure, indices_by_structure

def _slow_reader(index_map, data):
    """
    Loop over each valid pixel in the index_map and add its coordinates and
    data to the flux_by_structure and indices_by_structure dicts

    This is slower than _fast_reader but faster than that implementation would
    be without find_objects.  The bottleneck is doing `index_map == idx` N
    times.
    """
    flux_by_structure, indices_by_structure = {},{}
    # Do a fast iteration through d.data, adding the indices and data values
    # to the two dictionaries declared above:
    indices = np.array(np.where(index_map > -1)).transpose()

    log.debug('Creating index maps for {0} coordinates...'.format(len(indices)))
    for coord in ProgressBar(indices):
        coord = tuple(coord)
        idx = index_map[coord]
        if idx in flux_by_structure:
            flux_by_structure[idx].append(data[coord])
            indices_by_structure[idx].append(coord)
        else:
            flux_by_structure[idx] = [data[coord]]
            indices_by_structure[idx] = [coord]

    return flux_by_structure, indices_by_structure

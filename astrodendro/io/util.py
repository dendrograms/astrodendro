import numpy as np

from .. import six
from astropy.utils.console import ProgressBar
from astropy import log


def newick_from_json(d):
    """
    If "d" is a JSON-derived dict (e.g., `json.load('my_newick.json')`), this
    will convert all of the keys from strings to integers, resulting in
    something with identical format to `parse_newick`'s result.

    Example:
        >>> newick = ''.join(chr(x) for x in hdus[3].data.flat)
        >>> rslt = astrodendro.io.util.parse_newick(newick)
        >>> jj = json.dumps(rslt)
        >>> rdtrip = newick_from_json(json.loads(jj))
        >>> rdtrip == rslt
        True
    """
    new = {}
    for k, v in d.items():
        new_v = v
        if isinstance(v, dict):
            new_v = newick_from_json(v)
        elif isinstance(v, list):
            new_v = list()
            for x in v:
                if isinstance(x, dict):
                    new_v.append(newick_from_json(x))
                else:
                    new_v.append(x)
            new_v = tuple(new_v)
        new[int(k)] = new_v
    return new

def parse_newick(string):
    items = {}

    # Find maximum level
    current_level = 0
    max_level = 0
    log.debug("String starts with {0}".format(string[:100]))
    log.debug("String loading... newick has {0} chars, {1} ('s"
              .format(len(string), string.count("(")))
    for i, c in enumerate(string):
        if c == '(':
            current_level += 1
        if c == ')':
            current_level -= 1
        max_level = max(max_level, current_level)

    # Loop through levels and construct tree
    log.debug('Tree loading... max_level={0}'.format(max_level))
    for level in ProgressBar(range(max_level, 0, -1)):

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

    def collect(d):
        for item in d:
            if item in items:
                collect(items[item])
                d[item] = (items[item], d[item])
        return

    collect(items['trunk'])

    return items['trunk']


def _construct_tree(dend, repr, indices_by_structure, flux_by_structure):
    from ..structure import Structure
    structures = []
    for idx in repr:
        idx = int(idx)
        structure_indices = indices_by_structure[idx]
        f = flux_by_structure[idx]
        if type(repr[idx]) == tuple:
            sub_structures_repr = repr[idx][0]  # Parsed representation of sub structures
            sub_structures = _construct_tree(dend,
                                             sub_structures_repr,
                                             indices_by_structure,
                                             flux_by_structure)
            for i in sub_structures:
                dend._structures_dict[i.idx] = i
            b = Structure(structure_indices, f, children=sub_structures, idx=idx, dendrogram=dend)
            # Correct merge levels - complicated because of the
            # order in which we are building the tree.
            # What we do is look at the heights of this branch's
            # 1st child as stored in the newick representation, and then
            # work backwards to compute the merge level of this branch
            #
            # these five lines were not used
            #first_child_repr = six.next(six.itervalues(sub_structures_repr))
            #if type(first_child_repr) == tuple:
            #    height = first_child_repr[1]
            #else:
            #    height = first_child_repr
            dend._structures_dict[idx] = b
            structures.append(b)
        else:
            ell = Structure(structure_indices, f, idx=idx, dendrogram=dend)
            structures.append(ell)
            dend._structures_dict[idx] = ell
    return structures

def parse_dendrogram(parsed_newick, data, index_map, params, wcs=None):
    from ..dendrogram import Dendrogram

    d = Dendrogram()
    d.ndim = len(data.shape)

    d._structures_dict = {}
    d.data = data
    d.index_map = index_map
    d.params = params
    d.wcs = wcs

    try:
        flux_by_structure, indices_by_structure = _fast_reader(d.index_map, data)
    except ImportError:
        flux_by_structure, indices_by_structure = _slow_reader(d.index_map, data)

    log.debug('Constructing tree...')
    d.trunk = _construct_tree(d, parsed_newick, indices_by_structure, flux_by_structure)
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

    # find_objects returns a tuple that includes many None values that we
    # need to get rid of.
    object_slices = [x for x in object_slices if x is not None]

    index_cube = np.indices(index_map.shape)

    # Need to have same length, otherwise assumptions above are wrong
    assert len(idxs) == len(object_slices)
    log.debug('Creating index maps for {0} indices...'.format(len(idxs)))

    p = ProgressBar(len(object_slices))
    for idx,sl in zip(idxs, object_slices):
        match = index_map[sl] == idx
        sl2 = (slice(None),) + sl
        match_inds = index_cube[sl2][:, match]
        #coords = list(zip(*match_inds))
        dd = data[sl][match].tolist()
        flux_by_structure[idx] = dd
        indices_by_structure[idx] = match_inds.T
        p.update()

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

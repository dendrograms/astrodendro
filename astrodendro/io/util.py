import numpy as np

from .. import six


def parse_newick(string):

    items = {}

    # Find maximum level
    current_level = 0
    max_level = 0
    for i, c in enumerate(string):
        if c == '(':
            current_level += 1
        if c == ')':
            current_level -= 1
        max_level = max(max_level, current_level)

    # Loop through levels and construct tree
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

    # Do a fast iteration through d.data, adding the indices and data values
    # to the two dictionaries declared above:
    indices = np.indices(d.data.shape).reshape(d.data.ndim, np.prod(d.data.shape)).transpose()

    for coord in indices:
        coord = tuple(coord)
        idx = d.index_map[coord]
        if idx > -1:
            try:
                flux_by_structure[idx].append(d.data[coord])
                indices_by_structure[idx].append(coord)
            except KeyError:
                flux_by_structure[idx] = [d.data[coord]]
                indices_by_structure[idx] = [coord]

    d.trunk = _construct_tree(parse_newick(newick))
    # To make the structure.level property fast, we ensure all the items in the
    # trunk have their level cached as "0"
    for structure in d.trunk:
        structure._level = 0  # See the @property level() definition in structure.py

    d._index()
    return d

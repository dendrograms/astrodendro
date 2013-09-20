# Licensed under an MIT open source license - see LICENSE

# Notes:
# - A structure is a Leaf or a Branch
# - An ancestor is the largest structure that a structure is part of

import numpy as np

from .structure import Structure
from .progressbar import AnimatedProgressBar
from .io import IO_FORMATS
from . import pruning
from . import six


def _sorted_by_idx(d):
    return sorted(d, key=lambda s: s.idx)


class Dendrogram(object):
    """
    This class is used to compute and represent a dendrogram for a given
    dataset. To create a dendrogram from an array, use the
    :meth:`~astrodendro.dendrogram.Dendrogram.compute` class method::

        >>> from astrodendro import Dendrogram
        >>> d = Dendrogram.compute(array)

    Once the dendrogram has been computed, you can explore it programmatically
    using the ``trunk`` attribute, which allows you to access the base-level
    structures in the dendrogram::

        >>> d.trunk
        [<Structure type=leaf idx=101>,
         <Structure type=branch idx=2152>,
         <Structure type=leaf idx=733>,
         <Structure type=branch idx=303>]

    Structures can then be recursively explored. For more information on
    attributes and methods available for structures, see the
    :class:`~astrodendro.structure.Structure` class.

    The dendrogram can also be explored using an interactive viewer. To use
    this, use the :meth:`~astrodendro.dendrogram.Dendrogram.viewer` method::

        >>> d.viewer()

    and an interactive Matplotlib window should open.

    Finally, the :meth:`~astrodendro.dendrogram.Dendrogram.plotter` method can
    be used to facilitate the creation of plots:

        >>> p = d.plotter()

    For more information on using the plotter and other aspects of the
    :class:`~astrodendro.dendrogram.Dendrogram` class, see the online
    documentation.
    """

    def __init__(self):
        self.data = None
        self.n_dim = 0
        # Put in a friendly error message to make sure nobody confuses the
        # static methods for creating a dendrogram with instance methods:

        def static_warning(self, *args, **kwargs):
            err = "Invalid use of static method. Try d=Dendrogram.compute(data)"
            err += " or d=Dendrogram.load_from(file)"
            raise AttributeError(err)
        self.compute = static_warning
        self.load_from = static_warning

    @staticmethod
    def compute(data, min_value=-np.inf, min_delta=0, min_npix=0,
                is_independent=None, verbose=False):
        """
        Compute a dendrogram from a Numpy array.

        Parameters
        ----------
        data : `~numpy.ndarray`
            The n-dimensional array to compute the dendrogram for
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

        Examples
        --------

        The following example demonstrates how to compute a dendrogram from an
        dataset contained in a FITS file::

            >>> from astropy.io import fits
            >>> array = fits.getdata('observations.fits')
            >>> from astrodendro import Dendrogram
            >>> d = Dendrogram.compute(array)

        Notes
        -----
        More information about the above parameters is available from the
        online documentation at [www.dendrograms.org](www.dendrograms.org).
        """
        tests = [pruning.min_delta(min_delta),
                 pruning.min_npix(min_npix)]
        if is_independent is not None:
            if hasattr(is_independent, '__iter__'):
                tests.extend(is_independent)
            else:
                tests.append(is_independent)
        is_independent = pruning.all_true(tests)

        self = Dendrogram()
        self.data = data
        self.n_dim = len(data.shape)
        # For reference, store the parameters used:
        self.params = dict(min_npix=min_npix, min_value=min_value,
                           min_delta=min_delta)

        # Create a list of all points in the cube above min_value
        keep = self.data > min_value
        data_values = self.data[keep]
        indices = np.vstack(np.where(keep)).transpose()

        if verbose:
            print("Generating dendrogram using {:,} of {:,} pixels ({}% of data)".format(data_values.size, self.data.size, (100 * data_values.size / self.data.size)))
            progress_bar = AnimatedProgressBar(end=max(data_values.size, 1), width=40, fill='=', blank=' ')

        # Define index array indicating what structure each cell is part of
        # We expand each dimension by one, so the last value of each
        # index (accessed with e.g. [nx,#,#] or [-1,#,#]) is always zero
        # This permits an optimization below when finding adjacent structures
        self.index_map = np.zeros(np.add(self.data.shape, 1), dtype=np.int32)

        # Dictionary of currently-defined structures:
        structures = {}

        # Define a list of offsets we add to any coordinate to get the indices
        # of all neighbouring pixels
        if self.n_dim == 3:
            neighbour_offsets = np.array([(0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0)])
        elif self.n_dim == 2:
            neighbour_offsets = np.array([(0, -1), (0, 1), (-1, 0), (1, 0)])
        elif self.n_dim == 1:
            neighbour_offsets = np.array([(-1, ), (1, ), ])
        else:  # N-dimensional case. Analogous to the above.
            neighbour_offsets = np.concatenate((
                np.identity(self.n_dim, dtype=int),
                np.identity(self.n_dim, dtype=int) * -1
            ))

        # Loop from largest to smallest data_value value. Each time, check if
        # the pixel connects to any existing leaf. Otherwise, create new leaf.

        count = 0

        for i in np.argsort(data_values)[::-1]:

            def next_idx():
                return i + 1
                # Generate IDs index i. We add one to avoid ID 0

            data_value = data_values[i]
            coord = tuple(indices[i])

            # Print stats
            count += 1
            if verbose and (count % 100 == 0):
                progress_bar + 100
                progress_bar.show_progress()

            # Check if point is adjacent to any leaf
            # We don't worry about the edges, because overflow or underflow in
            # any one dimension will always land on an extra "padding" cell
            # with value zero added above when index_map was created
            indices_adjacent = [tuple(c) for c in np.add(neighbour_offsets, indices[i])]
            adjacent = [self.index_map[c] for c in indices_adjacent if self.index_map[c]]

            # Replace adjacent elements by its ancestor
            adjacent = [structures[a].ancestor for a in adjacent]

            # Remove duplicates
            adjacent = _sorted_by_idx(set(adjacent))

            # What happens next depends on how many unique adjacent structures there are

            if not adjacent:  # No adjacent structures;  Create new leaf:

                # Set absolute index of the new element
                idx = next_idx()

                # Create leaf
                leaf = Structure(coord, data_value, idx=idx, dendrogram=self)

                # Add leaf to overall list
                structures[idx] = leaf

                # Set absolute index of pixel in index map
                self.index_map[coord] = idx

            elif len(adjacent) == 1:  # Add to existing leaf or branch

                # Add point to structure
                adjacent[0]._add_pixel(coord, data_value)

                # Set absolute index of pixel in index map
                self.index_map[coord] = adjacent[0].idx

            else:  # Merge leaves

                # At this stage, the adjacent structures might consist of an
                # arbitrary number of leaves and branches.

                # Find all leaves that are not important enough to be
                # kept separate. These leaves will now be treated the
                # same as the pixel under consideration
                merge = [structure for structure in adjacent
                         if structure.is_leaf and
                         (structure.vmax == data_value or
                          not is_independent(structure, index=coord,
                                            value=data_value))]

                # Remove merges from list of adjacent structures
                for structure in merge:
                    adjacent.remove(structure)

                # Now, figure out what object this pixel belongs to
                # How many significant adjacent structures are left?

                if not adjacent:  # if len(adjacent) == 0:
                    # There are no separate leaves left (and no branches), so pick the
                    # first one as the reference and merge all the others onto it
                    belongs_to = merge.pop()
                    belongs_to._add_pixel(coord, data_value)
                elif len(adjacent) == 1:
                    # There is one significant adjacent leaf/branch left.
                    belongs_to = adjacent[0]
                    belongs_to._add_pixel(coord, data_value)
                else:
                    # Create a branch
                    belongs_to = Structure(coord, data_value,
                                           children=adjacent, idx=next_idx(),
                                           dendrogram=self)
                    # Add branch to overall list
                    structures[belongs_to.idx] = belongs_to

                # Set absolute index of pixel in index map
                self.index_map[coord] = belongs_to.idx

                # Add all insignificant leaves in 'merge' to the same object as this pixel:
                for m in merge:
                    # print "Merging leaf %i onto leaf %i" % (i, idx)
                    # Remove leaf
                    structures.pop(m.idx)
                    # Merge the insignificant structure that this pixel now belongs to:
                    belongs_to._merge(m)
                    # Update index map
                    m._fill_footprint(self.index_map, belongs_to.idx)

        if verbose:
            progress_bar.progress = 100  # Done
            progress_bar.show_progress()
            print("")  # newline

        # Create trunk from objects with no ancestors
        self.trunk = _sorted_by_idx([structure for structure in six.itervalues(structures) if structure.parent is None])

        # Remove orphan leaves that aren't large enough
        leaves_in_trunk = [structure for structure in self.trunk if structure.is_leaf]
        for leaf in leaves_in_trunk:
            if not is_independent(leaf):
                # This leaf is an orphan, so remove all references to it:
                structures.pop(leaf.idx)
                self.trunk.remove(leaf)
                leaf._fill_footprint(self.index_map, 0)

        # To make the structure.level property fast, we ensure all the structures in the
        # trunk have their level cached as "0"
        for structure in self.trunk:
            structure._level = 0  # See the definition of level() in structure.py

        # Save a list of all structures accessible by ID
        self._structures_dict = structures

        #remove border from index map
        s = tuple(slice(0, s, 1) for s in data.shape)
        self.index_map = self.index_map[s]

        self._index()

        # Return the newly-created dendrogram:
        return self


    def _index(self):
        # add dendrogram index
        ti = TreeIndex(self)

        for s in six.itervalues(self._structures_dict):
            s._tree_index = ti


    @property
    def trunk(self):
        """
        A list of all structures that have no parent structure and form the
        base of the tree.
        """
        return self._trunk

    @trunk.setter
    def trunk(self, value):
        self._trunk = value

    @staticmethod
    def load_from(filename, format=None):
        """
        Load a previously computed dendrogram from a file.

        Parameters
        ----------
        filename : str
            The name of the file to load the dendrogram from. By default, the
            file format will be automatically detected from the file
            extension. At this time, only HDF5 files (extension ``.hdf5``) are
            supported.
        format : str, optional
            The format to use to read the file. By default, this is not used
            and the format is auto-detected from the file extension. At this
            time, the only format supported is ``'hdf5'``.
        """
        from .io import load_dendrogram
        return load_dendrogram(filename, format=format)

    def save_to(self, filename, format=None):
        """
        Save the dendrogram to a file.

        Parameters
        ----------
        filename : str
            The name of the file to save the dendrogram to. By default, the
            file format will be automatically detected from the file
            extension. At this time, only HDF5 files (extension ``.hdf5``) are
            supported.
        format : str, optional
            The format to use for the file. By default, this is not used and
            the format is auto-detected from the file extension. At this time,
            the only format supported is ``'hdf5'``.
        """
        from .io import save_dendrogram
        return save_dendrogram(self, filename, format=format)

    @property
    def leaves(self):
        """
        A flattened list of all leaves in the dendrogram
        """
        return [i for i in six.itervalues(self._structures_dict) if i.is_leaf]

    def to_newick(self):
        #this caches newicks, and prevents too much recursion
        [s.newick for s in reversed(list(self.all_structures))]

        return "(%s);" % ','.join([structure.newick for structure
                                   in self.trunk])

    def structure_at(self, indices):
        """
        Get the structure at the specified pixel coordinate.

        This will return None if no structure includes the specified pixel
        coordinates.
        """
        idx = self.index_map[indices]
        if idx:
            return self._structures_dict[idx]
        return None

    @property
    def all_structures(self):
        """
        Yields an iterator over all structures in the dendrogram, in prefix order.
        """

        todo = list(self.trunk)
        while len(todo) > 0:
            st = todo.pop(0)
            yield st
            todo = st.children + todo

    def __getitem__(self, key):
        """Fetch structures by index value"""
        return self._structures_dict[key]

    def __len__(self):
        """Return number of structures in dendrogram"""
        return len(self._structures_dict)

    def __iter__(self):
        return self.all_structures

    def __eq__(self, other):
        if not isinstance(other, Dendrogram):
            return False

        if not (self.data == other.data).all():
            return False

        # structures should have the same extent,
        # but idx values need not be identical. This
        # tests the index map for that
        u, ind = np.unique(self.index_map, return_index=True)
        u, ind2 = np.unique(self.index_map, return_index=True)
        return (np.sort(ind) == np.sort(ind2)).all()

    def plotter(self):
        """
        Return a :class:`~astrodendro.plot.DendrogramPlotter` instance that makes it easier to construct plots.
        """
        from .plot import DendrogramPlotter
        return DendrogramPlotter(self)

    def viewer(self):
        """
        Launch an interactive viewer to explore the dendrogram.

        This functionality is only available for 2- or 3-d datasets.
        """
        from .viewer import BasicDendrogramViewer
        return BasicDendrogramViewer(self)


class TreeIndex(object):
    def __init__(self, dendrogram):
        """
        Object that efficiently extracts
        the locations of Structures in an ndarray

        Parameters
        ----------
        dendrogram : Dendrogram instance
                    The dendrogram to index
        """
        index_map = dendrogram.index_map
        trunk = dendrogram.trunk

        sz = index_map.size
        nd = len(index_map.shape)

        assert sz == dendrogram.data.size
        assert index_map.min() >= 0

        #map ids to [0, 1, ...] for storage efficiency
        uniq, bins = np.unique(index_map, return_inverse=True)
        packed = dict((u, i) for i, u in enumerate(uniq))

        flat_idx = index_map.ravel()
        ri = np.argsort(bins)
        idx_ct = np.bincount(bins)
        idx_sub_ct = {}
        idx_cdf = np.hstack((0, np.cumsum(idx_ct)))

        #efficiently build up npix values
        structures = reversed(sorted(dendrogram._structures_dict.values(),
                                key=lambda x: x.level))
        for st in structures:
            idx_sub_ct[st.idx] = idx_ct[packed[st.idx]]
            idx_sub_ct[st.idx] += sum(idx_sub_ct[c.idx] for c in st.children)

        #build a 1D index array with the following properties
        # - values in index reference locations in flattened index_map
        # - every structure (+ subtree) is a continuous slice of index
        # - index[offset[i]] is the first location for (packed) structure i
        # - npix[i] gives number of pixels in (packed) structure i,
        #   exluding subtree
        # - npix_subtree[i] is like above, but includes subtrees
        #
        # In summary, the locations in the flattened_index map
        # for structure i excluding subrees is
        #    pi = packed[i]
        #    index[offset[pi] : offset[pi] + npix[pi]]
        # and including subtrees is
        #    index[offset[pi] : offset[pi] + npix_subtree[pi]]
        offset = np.zeros(idx_ct.size, dtype=np.int)
        npix = offset * 0
        npix_subtree = offset * 0

        index = np.zeros(sz, dtype=np.int)
        order = dendrogram.all_structures

        pos = 0
        for o in order:
            sid = packed[o.idx]
            offset[sid] = pos
            npix[sid] = idx_ct[sid]
            npix_subtree[sid] = idx_sub_ct[o.idx]
            idx = ri[idx_cdf[sid]: idx_cdf[sid] + npix[sid]]
            assert (flat_idx[idx] == o.idx).all()
            index[pos: pos + npix[sid]] = idx
            pos += npix[sid]

        #turn inds back into an ndim index
        self._index = tuple(n.ravel()[index] for n in
                            np.indices(index_map.shape))

        self._data = dendrogram.data
        self._offset = offset
        self._npix = npix
        self._npix_subtree = npix_subtree
        self.packed = packed

    def indices(self, sid, subtree=True):
        """
        Return pixel indices associated with a dendrogram structure

        Returns the pixels in the original dendrogram array
        which are associated with a particular structure id.

        Parameters
        ----------
        sid : integer
              The structure index to lookup. Stored in `Structure.idx`
        subtree : bool, optional
              If true, return indices for subtrees as well. Default=False

        Returns
        -------
        A tuple of integer ndarrays, akin to np.where().
        This can be directly used as an index into the dendrogram
        data array
        """
        sid = self.packed[sid]
        i0 = self._offset[sid]
        di = self._npix_subtree[sid] if subtree else self._npix[sid]
        return tuple(ind[i0: i0 + di] for ind in self._index)

    def values(self, sid, subtree=True):
        return self._data[self.indices(sid, subtree=subtree)]

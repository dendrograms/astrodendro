# Licensed under an MIT open source license - see LICENSE

import numpy as np


def cached_property(func):
    memo = {}
    def result(self):
        if self not in memo:
            memo[self] = func(self)
        return memo[self]

    return property(result)


def prefix_visit(s, key=None, reverse=False):
    todo = [s]
    while todo:
        st = todo.pop(0)
        yield st
        children = st.children
        if key is not None:
            children = sorted(children, key=key, reverse=reverse)
        todo = children + todo


class Structure(object):
    """
    A structure in the dendrogram, for example a leaf or a branch.

    A structure that is part of a dendrogram knows which other structures it is
    related to. For example, it is possible to get the parent structure
    containing the present structure ``s`` by using the ``parent`` attribute::

        >>> s.parent
        <Structure type=branch idx=2152>

    Likewise, the ``children`` attribute can be used to get a list of all
    sub-structures::

        >>> s.children
        [<Structure type=branch idx=1680>, <Structure type=branch idx=5771>]

    A number of attributes and methods are available to explore the structure
    in more detail, such as the ``indices`` and ``values`` methods, which
    return the indices and values of the pixels that are part of the
    structure. These and other methods have a ``subtree=`` option, which if
    ``True`` (the default) returns the quantities related to structure and all
    sub-structures, and if ``False`` includes only the pixels that are part of
    the structure, but excluding any sub-structure.
    """

    ###########################################################################
    #   The following methods are used only by the Dendrogram class,          #
    #   for computing the dendrogram and should never be used manually:       #
    ###########################################################################

    def __init__(self, indices, values, children=[], idx=None, dendrogram=None):

        self._dendrogram = dendrogram
        self.parent = None
        self.children = children

        # Make sure that the children have a reference to the present structure
        for child in children:
            child.parent = self

        if np.isscalar(values):  # values are for a single pixel
            self._indices = [indices]
            self._values = [values]
            self._vmin, self._vmax = values, values
        elif isinstance(indices, list) and isinstance(values, list):
            self._indices = indices
            self._values = values
            self._vmin, self._vmax = min(values), max(values)
        else:  # could be an array or iterator
            self._indices = [x for x in indices]
            self._values = [x for x in values]
            self._vmin, self._vmax = min(values), max(values)

        self._smallest_index = min(self._indices)

        self.idx = idx

        self._reset_cache()

    def _reset_cache(self):
        self._level = None
        self._ancestor = None
        self._descendants = None
        self._npix_total = None
        self._peak = None
        self._peak_subtree = None
        self._tree_index = None

    @property
    def smallest_index(self):
        return self._smallest_index

    @smallest_index.setter
    def smallest_index(self, value):
        self._smallest_index = value

    @property
    def parent(self):
        """
        The parent structure containing the present structure.
        """
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value

    @property
    def children(self):
        """
        A list of all the sub-structures contained in the present structure.
        """
        return self._children

    @children.setter
    def children(self, value):
        self._children = value

    @property
    def is_leaf(self):
        """
        Whether the present structure is a leaf.
        """
        return not self.children

    @property
    def is_branch(self):
        """
        Whether the present structure is a branch.
        """
        return not self.is_leaf

    def indices(self, subtree=True):
        """
        The indices of the pixels in this branch.

        Parameters
        ----------
        subtree : bool, optional
            Whether to recursively include all sub-structures
        """

        if self._tree_index is not None:
            return self._tree_index.indices(self.idx, subtree=subtree)

        if subtree:
            sub_indices = [self.indices(subtree=False)]
            for child in self.children:
                sub_indices.append(child.indices(subtree=True))
            return tuple(np.hstack(arrs) for arrs in zip(*(sub_indices)))
        else:
            return tuple(np.atleast_1d(i) for i in zip(*self._indices))

    def values(self, subtree=True):
        """
        The values of the pixels in this branch.

        Parameters
        ----------
        subtree : bool, optional
            Whether to recursively include all sub-structures
        """

        if self._tree_index is not None:
            return self._tree_index.values(self.idx, subtree=subtree)

        if subtree:
            sub_values = [self.values(subtree=False)]
            for child in self.children:
                sub_values.append(child.values(subtree=True))
            return np.hstack(sub_values)
        else:
            return np.atleast_1d(self._values)

    @property
    def vmin(self):
        """
        The minimum value of pixels belonging to the branch (excluding sub-structure).
        """
        return self._vmin

    @property
    def vmax(self):
        """
        The maximum value of pixels belonging to the branch (excluding sub-structure).
        """
        return self._vmax

    def _add_pixel(self, index, value):
        """
        Add a pixel to this leaf.

        This is a private method only intended for use by
        :meth:`~astrodendro.dendrogram.Dendrogram.compute`

        Parameters
        ----------
        index : tuple
            The pixel coordinates
        value : float
            The value of the pixel
        """
        self._indices.append(index)
        self._values.append(value)
        self._vmin, self._vmax = min(value, self.vmin), max(value, self.vmax)
        self._smallest_index = min(self._smallest_index, index)
        self._reset_cache()

    def _merge(self, structure):
        """
        Merge a structure with the present structure.

        This is a private method only intended for use by
        :meth:`~astrodendro.dendrogram.Dendrogram.compute`
        """
        self._indices.extend(structure._indices)
        self._values.extend(structure._values)
        self._vmin, self._vmax = min(structure.vmin, self.vmin), max(structure.vmax, self.vmax)
        self._smallest_index = min(structure._smallest_index, self._smallest_index)
        self._reset_cache()

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################

    @property
    def height(self):
        """
        This is defined as the minimum value in the children structures, or the
        peak value of the present structure if it has no children.
        """
        return min(c.vmin for c in self.children) if self.children else self.vmax

    @property
    def ancestor(self):
        """
        Find the ancestor of this leaf/branch non-recursively.

        This always returns an object (may return self if the object has no
        parent). Results are partially cached to reduce recursion depth. The
        caching assumes that once an object has been given a parent, that
        parent will never change.

        This method should be equivalent to:

            if self.parent == None:
                return self
            else:
                return self.parent.ancestor
        """

        if self.parent is None:
            return self

        if not self._ancestor:
            self._ancestor = self.parent

        # Use a loop rather than recursion to update the cached ancestor:
        while self._ancestor.parent:
            a = self._ancestor
            if a._ancestor:
                self._ancestor = a._ancestor
            else:
                self._ancestor = a.parent

        return self._ancestor

    ###########################################################################
    #   The following methods are only reliable after the entire tree is      #
    #   computed. They should not be used in dendrogram.py                    #
    ###########################################################################

    def _fill_footprint(self, array, level, recursive=True):
        """
        Set all corresponding points in `array` to the level of the structure.

        Parameters
        ----------
        array : :class:`~numpy.ndarray`
            The array to write the footprint to
        recursive : bool
            Whether to also add the footprint for the sub-structures
        """
        if recursive:
            for child in self.children:
                child._fill_footprint(array, level + 1)
        array[self.indices(subtree=False)] = level

    @property
    def level(self):
        """
        The level of the structure, i.e. how many structures need to be traversed to reach the present structure.

        This is 0 for structures in the trunk, with values increasing in steps of 1
        towards the leaves.
        """

        if self._level is None:
            if not self.parent:
                self._level = 0
            elif self.parent._level is not None:
                self._level = self.parent._level + 1
            else:
                # We could just use:
                #  self._level = self.parent.level + 1
                # But to avoid recursion, and keep things fast, we do it this way instead:
                obj = self.parent
                diff = 1
                while obj._level is None:
                    obj = obj.parent
                    diff += 1
                    # Note: we are counting on the dendrogram computation to
                    # ensure that ._level=0 for all structures in the trunk
                self._level = obj._level + diff
                self.parent._level = self._level - 1

        return self._level

    @cached_property
    def newick(self):
        """
        Newick representation of this structure.
        """
        if self.idx is None:
            raise ValueError("Cannot return Newick representation if idx is not set")
        if self.children:
            newick_items = [child.newick for child in self.children]
            return "(%s)%s:%.3f" % (','.join(newick_items), self.idx, self.height)
        else:
            return "%i:%.3f" % (self.idx, self.height)

    @property
    def descendants(self):
        """
        Get a flattened list of all child leaves and branches.
        """

        if self._descendants is None:
            self._descendants = []
            to_add = [self]  # branches with children we will need to add to the list
            while True:
                children = []
                list(map(children.extend, [branch.children for branch in to_add]))
                self._descendants.extend(children)
                # Then proceed, essentially recursing through child branches:
                to_add = [b for b in children if not b.is_leaf]
                if not to_add:
                    break

        return self._descendants

    def get_npix(self, subtree=True):
        """
        Return the number of pixels in this structure.

        Parameters
        ----------
        subtree : bool, optional
            Whether to recursively include all sub-structures when counting the
            pixels.

        Returns
        -------
        n_pix : int
            The number of pixels in this structure
        """

        if subtree:
            if self._npix_total is None:
                self._npix_total = self.values(subtree=True).size
            return self._npix_total
        else:
            return len(self.values(subtree=False))

    def get_peak(self, subtree=True):
        """
        Return (index, value) for the pixel with maximum value.

        Parameters
        ----------
        subtree : bool, optional
            Whether to recursively include all sub-structures when searching
            for the peak.

        Returns
        -------
        index : tuple
            The n-dimensional index of the peak pixel
        value : float
            The value of the peak pixel
        """
        # populate caches without recursion
        def key(x):
            return x[1]

        if self._peak is None:
            for s in reversed(list(prefix_visit(self))):
                s._peak = (s._indices[s._values.index(s.vmax)],
                          s.vmax)
                if s.is_leaf:
                    s._peak_subtree = s._peak
                else:
                    s._peak_subtree = max((c._peak_subtree
                                           for c in s.children),
                                           key=key)
                    s._peak_subtree = max(s._peak_subtree, s._peak, key=key)

        if not subtree:
            return self._peak

        return self._peak_subtree

    def __repr__(self):
        if self.is_leaf:
            return "<Structure type=leaf idx={0}>".format(self.idx)
        else:
            return "<Structure type=branch idx={0}>".format(self.idx)

    def sorted_leaves(self, sort_key=lambda s: s.get_peak(subtree=True)[1],
                          reverse=False, subtree=True):
        """
        Return a list of sorted leaves.

        Parameters
        ----------
        sort_key : function, optional
            A function which given a structure will return a scalar that is
            then used for sorting. By default, this is set to a function that
            returns the peak value of a structure (including descendants).
        reverse : bool, optional
            Whether to reverse the sorting.
        subtree : bool, optional
            Whether to recursively include all sub-structures in the list.

        Returns
        -------
        leaves : list
            A list of sorted leaves
        """
        if self.is_leaf:
            return [self]
        if not subtree:
            return [s for s in sorted(self.children, key=sort_key,
                                      reverse=reverse)
                    if s.is_leaf]
        # flip reverse keyword, so that nodes are properly
        # sorted after reversed()
        sts = reversed(list(prefix_visit(self, key=sort_key,
                                         reverse=not reverse)))
        return [s for s in sts if s.is_leaf]

    def get_mask(self, shape=None, subtree=True):
        """
        Return a boolean mask outlining the structure.

        Parameters
        ----------
        shape : tuple, optional
            The shape of the array upon which to compute the mask. This is only
            required if the structure is not attached to a dendrogram.
        subtree : bool, optional
            Whether to recursively include all sub-structures in the mask.

        Returns
        -------
        mask : :class:`~numpy.ndarray`
            The mask outlining the structure (``False`` values are used outside
            the structure, and ``True`` values inside).
        """
        if shape is None:
            shape = self._dendrogram.data.shape
        indices = self.indices(subtree=True) if subtree else self.indices
        mask = np.zeros(shape, dtype=bool)
        mask[indices] = True
        return mask

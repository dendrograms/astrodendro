# Computing Astronomical Dendrograms
# Copyright (c) 2011-2012 Thomas P. Robitaille and Braden MacDonald
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import numpy as np


class Structure(object):
    """
    A structure in the dendrogram, for example a leaf or a branch.
    """

    ###########################################################################
    #   The following methods are used only by the Dendrogram class,          #
    #   for computing the dendrogram and should never be used manually:       #
    ###########################################################################

    def __init__(self, indices, values, children=[], idx=None):

        self.parent = None
        self.children = children

        # Make sure that the children have a reference to the present structure
        for child in children:
            child.parent = self

        if np.isscalar(values):  # values are for a single pixel
            self._indices = [indices]
            self._values = [values]
            self._vmin, self._vmax = values, values
        else:  # values are for a sequence of pixels
            self._indices = indices
            self._values = values
            self._vmin, self._vmax = min(values), max(values)

        self.idx = idx

        self._reset_cache()

    def _reset_cache(self):
        self._level = None
        self._ancestor = None
        self._descendants = None
        self._npix_total = None
        self._peak = None
        self._tree_index = None

    @property
    def is_leaf(self):
        """
        Whether the present structure is a leaf
        """
        return not self.children

    @property
    def is_branch(self):
        """
        Whether the present structure is a branch
        """
        return not self.is_leaf

    @property
    def indices_all(self):
        """
        The indices of the pixels in the branch, and sub-structures
        """
        # We only need to look at children, not all children structures,
        # since child.indices will include all children recursively.
        if self._tree_index is not None:
            return self._tree_index.indices(self.idx, subtree=True)

        sub_indices = [self.indices]
        for child in self.children:
            sub_indices.append(child.indices_all)
        return tuple(np.hstack(arrs) for arrs in zip(*(sub_indices)))

    @property
    def indices(self):
        """
        The indices of the pixels in this branch, excluding sub-structures
        """
        if self._tree_index is not None:
            return self._tree_index.indices(self.idx, subtree=False)

        return tuple(np.atleast_1d(i) for i in zip(*self._indices))

    @property
    def values_all(self):
        """
        The values of the pixels in the branch, and sub-structures
        """
        # We only need to look at children, not all children structures,
        # ``values`` for all children will include all children recursively.
        if self._tree_index is not None:
            return self._tree_index.values(self.idx, subtree=True)

        sub_values = [self.values]
        for child in self.children:
            sub_values.append(child.values_all)
        return np.hstack(sub_values)

    @property
    def values(self):
        """
        The values of the pixels in this branch, excluding sub-structures
        """
        if self._tree_index is not None:
            return self._tree_index.values(self.idx, subtree=False)

        return np.atleast_1d(self._values)

    @property
    def vmin(self):
        """
        The minimum value of pixels belonging to the branch (excluding sub-structure)
        """
        return self._vmin

    @property
    def vmax(self):
        """
        The maximum value of pixels belonging to the branch (excluding sub-structure)
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
        self._reset_cache()

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################

    @property
    def height(self):
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

    def fill_footprint(self, array, level, recursive=True):
        """
        Set all corresponding points in `array` to the level of the structure

        Parameters
        ----------
        array : `~numpy.ndarray`
            The array to write the footprint to
        recursive : bool
            Whether to also add the footprint for the sub-structures
        """
        if recursive:
            for child in self.children:
                child.fill_footprint(array, level + 1)
        array[self.indices] = level

    @property
    def level(self):
        """
        The level of the structure.

        This is 0 for nodes in the trunk, with values increasing in steps of 1
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
                    # ensure that ._level=0 for all nodes in the trunk
                self._level = obj._level + diff
                self.parent._level = self._level - 1

        return self._level

    @property
    def newick(self):
        """
        Newick representation of this structure
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
                map(children.extend, [branch.children for branch in to_add])
                self._descendants.extend(children)
                # Then proceed, essentially recursing through child branches:
                to_add = [b for b in children if not b.is_leaf]
                if not to_add:
                    break

        return self._descendants

    def get_npix(self, subtree=False):
        """
        Return the number of pixels in this structure.

        If subtree is True, the result is a sum that includes all child nodes.
        """

        if not subtree:
            return len(self.values)
        else:
            if self._npix_total is None:
                self._npix_total = self.values_all.size
            return self._npix_total

    def get_peak(self, subtree=False):
        """
        Return (indices, values) for the pixel with maximum value

        If subtree=True, will search all descendant nodes too.
        """
        if self._peak is None:
            self._peak = (self._indices[self._values.index(self.vmax)], self.vmax)
            # Note the above cached value never includes descendants

        if not subtree:
            return self._peak
        else:
            found = self._peak
            for node in self.descendants:
                if found[1] < node.vmax:
                    found = node.get_peak()
            return found

    def __repr__(self):
        if self.is_leaf:
            return "<Structure type=leaf idx={0}>".format(self.idx)
        else:
            return "<Structure type=branch idx={0}>".format(self.idx)

    def get_sorted_leaves(self, sort_key=lambda s: s.get_peak(subtree=True)[1], reverse=False, subtree=False):
        if self.is_leaf:
            return [self]
        leaves = []
        for structure in sorted(self.children, key=sort_key, reverse=reverse):
            if structure.is_leaf:
                leaves.append(structure)
            elif subtree:
                leaves += structure.get_sorted_leaves(sort_key=sort_key, reverse=reverse, subtree=subtree)
        return leaves

    def get_mask(self, shape, subtree=False):
        indices = self.indices_all if subtree else self.indices
        mask = np.zeros(shape, dtype=bool)
        mask[indices] = True
        return mask

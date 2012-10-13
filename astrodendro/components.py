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

class Leaf(object):

    ###########################################################################
    #   The following methods are used only by the Dendrogram class,          #
    #   for computing the dendrogram and should never be used manually:       #
    ###########################################################################
    def __init__(self, coord, f, idx=None):
        if not hasattr(f, '__iter__'):
            # Normal initialization - coord and f are for a single pixel:
            self.coords = [coord]
            self.f = [f]
            self.fmin, self.fmax = f, f
        else:
            self.coords = coord
            self.f = f
            self.fmin, self.fmax = min(f), max(f)
        self.idx = idx
        self.parent = None
        self._ancestor = None # Cached ancestor, if any
        self._level = None # Cached "level" property - see below

    def _add_pixel(self, coord, f):
        """
        Add the pixel at `coord` with intensity `f` to this leaf
        This is only intended for use by Dendrogram.compute()
        """
        self.coords.append(coord)
        self.f.append(f)
        self.fmin, self.fmax = min(f, self.fmin), max(f, self.fmax)

    def _merge(self, leaf):
        """
        Pixels from the given leaf are now to be part of this leaf.
        This is only intended for use by Dendrogram.compute()
        """
        self.coords.extend(leaf.coords)
        self.f.extend(leaf.f)
        self.fmin, self.fmax = min(leaf.fmin, self.fmin), max(leaf.fmax, self.fmax)

    def add_footprint(self, image, level):
        "Fill in a map which shows the depth of the tree"
        for c in self.coords:
            image[c] = level

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################

    @property
    def height(self):
        if self.parent == None:
            return self.fmax - self.fmin
        else:
            return self.fmax - self.parent.merge_level
    
    @property
    def ancestor(self):
        """
        Find the ancestor of this leaf/branch non-recursively. Always returns
        an object (may return self if the object has no parent).
        Results are partially cached to reduce recursion depth.
        The caching assumes that once an object has been given a parent, that
        parent will never change.
        
        This method should be equivalent to:
            if self.parent == None:
                return self
            else:
                return self.parent.ancestor
        """
        if self.parent == None:
            return self
        if not self._ancestor:
            self._ancestor = self.parent
        while self._ancestor.parent:
            # Use a loop rather than recursion to update
            # the cached ancestor, if needed:
            a = self._ancestor
            if a._ancestor:
                self._ancestor = a._ancestor # Update our cached value
            else:
                self._ancestor = a.parent
        return self._ancestor
    
    ###########################################################################
    #   The following methods are only reliable after the entire tree is      #
    #   computed. They should not be used in dendrogram.py                    #
    ###########################################################################

    @property
    def level(self):
        " Level: 0 for nodes in the trunk, 1 for their immediate children, etc"
        if self._level is None:
            if not self.parent:
                self._level = 0
            elif self.parent._level is not None:
                self._level = self.parent._level + 1
            else:
                #We could just use:
                #  self._level = self.parent.level + 1
                #But to avoid recursion, and keep things fast, we do it this way instead:
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
        " Newick representation of this Leaf " 
        return "%i:%.3f" % (self.idx, self.height)

    def get_npix(self, descend=False):
        """
        Return the number of pixels in this Leaf.
        'descend' is ignored and is only for compatibility with Branch
        """
        return len(self.f)

    def get_peak(self, descend=False):
        """
        Return (coordinates, intensity) for the pixel with maximum value.
        'descend' is ignored and is only for compatibility with Branch
        """
        if not hasattr(self, '_peak'):
            self._peak = (self.coords[self.f.index(self.fmax)], self.fmax) 
        return self._peak

class Branch(Leaf):

    def __init__(self, children, coord, f, idx=None):
        self.merge_level = f # Record the exact flux level that triggered creation of this branch
        self.children = children
        for child in children:
            child.parent = self
        Leaf.__init__(self, coord, f, idx=idx)
        self._descendants = None # Cached value is not initially set

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################
    
    def add_footprint(self, image, level, recursive=True):
        if recursive:
            for child in self.children:
                child.add_footprint(image, level + 1)
        Leaf.add_footprint(self, image, level)
    
    
    ###########################################################################
    #   The following methods are only reliable after the entire tree is      #
    #   computed. They should not be used in dendrogram.py                    #
    ###########################################################################

    @property
    def newick(self):
        newick_items = [child.newick for child in self.children]
        return "(%s)%s:%.3f" % (','.join(newick_items), self.idx, self.height)
    
    @property
    def descendants(self):
        "Get a flattened list of all child leaves and branches. Non-recursive."
        if self._descendants is None:
            self._descendants = []
            to_add = [self] # branches with children we will need to add to the list
            while True:
                children = []
                map(children.extend, [branch.children for branch in to_add])
                self._descendants.extend(children)
                # Then proceed, essentially recursing through child branches:
                to_add = [b for b in children if type(b) is Branch]
                if not to_add:
                    break
        return self._descendants
    
    def get_npix(self, descend=False):
        """
        Return the number of pixels in this Branch.
        If descend=True, the result is a sum that includes all child nodes. 
        """
        if not descend:
            return len(self.f)
        # descend is True, so return total npix including all child nodes:
        if not hasattr(self, '_npix_total'): # _npix_total has not been cached
            self._npix_total = len(self.f)
            self._npix_total += sum([len(node.f) for node in self.descendants])
        return self._npix_total
    
    def get_peak(self, descend=False):
        """
        Return (coordinates, intensity) for the pixel with maximum value
        If descend=True, will search all descendant nodes too. 
        """
        if not hasattr(self, '_peak'):
            self._peak = (self.coords[self.f.index(self.fmax)], self.fmax)
            # Note the above cached value never includes descendants
        if not descend:
            return self._peak
        else:
            found = self._peak
            for node in self.descendants:
                if found[1] < node.fmax:
                    found = node.get_peak()
            return found
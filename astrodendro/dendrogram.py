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

# Notes:
# - A node is a Leaf or a Branch
# - An ancestor is the largest structure that a node is part of

import numpy as np
from .components import Branch, Leaf
from .newick import parse_newick
from .progressbar import AnimatedProgressBar


class Dendrogram(object):

    def __init__(self):
        self.data = None
        self.n_dim = 0

    @staticmethod
    def compute(data, min_intensity=-np.inf, min_npix=0, min_delta=0, verbose=False):
        self = Dendrogram()
        self.data = data
        self.n_dim = len(data.shape)
        # For reference, store the parameters used:
        self.min_intensity, self.min_npix, self.min_delta = min_intensity, min_npix, min_delta
        
        # Create a list of all points in the cube above min_intensity
        keep = self.data.ravel() > min_intensity
        intensity_values = self.data.ravel()[keep]
        coords = np.array(np.unravel_index( np.arange(self.data.size)[keep] , self.data.shape)).transpose()
        
        if verbose:
            print("Generating dendrogram using {:,} of {:,} pixels ({}% of data)".format(intensity_values.size, self.data.size, (100*intensity_values.size/self.data.size)))
            progress_bar = AnimatedProgressBar(end=max(intensity_values.size,1), width=40, fill='=', blank=' ')
            
        # Define index array indicating what node each cell is part of
        # We expand each dimension by one, so the last value of each
        # index (accessed with e.g. [nx,#,#] or [-1,#,#]) is always zero
        # This permits an optimization below when finding adjacent nodes
        self.index_map = np.zeros(np.add(self.data.shape, 1), dtype=np.int32)

        # Dictionary of currently-defined nodes:
        nodes = {}
        
        # Define a list of offsets we add to any coordinate to get the coords
        # of all neighbouring pixels
        if self.n_dim == 3:
            neighbour_offsets = np.array([(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)])
        elif self.n_dim == 2:
            neighbour_offsets = np.array([(0,-1),(0,1),(-1,0),(1,0)])
        elif self.n_dim == 1:
            neighbour_offsets = np.array([(-1,),(1,),])
        else: # N-dimensional case. Analogous to the above.
            neighbour_offsets = np.concatenate((
                np.identity(self.n_dim, dtype=int),
                np.identity(self.n_dim, dtype=int)*-1
            ))

        # Loop from largest to smallest intensity value. Each time, check if
        # the pixel connects to any existing leaf. Otherwise, create new leaf.
        
        count = 0

        for i in np.argsort(intensity_values)[::-1]:
            
            def next_idx():
                return i+1
                # Generate IDs index i. We add one to avoid ID 0
            
            intensity = intensity_values[i]
            coord = tuple(coords[i])
            
            # Print stats
            count += 1
            if verbose and (count % 100 == 0):
                progress_bar + 100
                progress_bar.show_progress()

            # Check if point is adjacent to any leaf
            # We don't worry about the edges, because overflow or underflow in 
            # any one dimension will always land on an extra "padding" cell 
            # with value zero added above when index_map was created
            indices_adjacent = [tuple(c) for c in np.add(neighbour_offsets, coords[i])]
            adjacent = [self.index_map[c] for c in indices_adjacent if self.index_map[c]]
            
            # Replace adjacent elements by its ancestor
            adjacent = [ nodes[a].ancestor for a in adjacent]

            # Remove duplicates
            adjacent = list(set(adjacent))

            # What happens next depends on how many unique adjacent structures there are

            if not adjacent:  # No adjacent structures;  Create new leaf:

                # Set absolute index of the new element
                idx = next_idx()

                # Create leaf
                leaf = Leaf(coord, intensity, idx=idx)

                # Add leaf to overall list
                nodes[idx] = leaf

                # Set absolute index of pixel in index map
                self.index_map[coord] = idx

            elif len(adjacent) == 1:  # Add to existing leaf or branch

                # Add point to node
                adjacent[0]._add_pixel(coord, intensity)

                # Set absolute index of pixel in index map
                self.index_map[coord] = adjacent[0].idx

            else:  # Merge leaves

                # At this stage, the adjacent nodes might consist of an 
                # arbitrary number of leaves and branches.

                # Find all leaves that are not important enough to be kept
                # separate. These leaves will now be treated the same as the pixel
                # under consideration
                merge = [node for node in adjacent
                         if type(node) is Leaf and
                         (node.fmax - intensity < min_delta or 
                          len(node.f) < min_npix or node.fmax == intensity)]

                # Remove merges from list of adjacent nodes
                for node in merge:
                    adjacent.remove(node)

                # Now, figure out what object this pixel belongs to
                # How many significant adjacent nodes are left?

                if not adjacent: #if len(adjacent) == 0:
                    # There are no separate leaves left (and no branches), so pick the
                    # first one as the reference and merge all the others onto it
                    belongs_to = merge.pop()
                    belongs_to._add_pixel(coord, intensity)
                elif len(adjacent) == 1:
                    # There is one significant adjacent leaf/branch left.
                    belongs_to = adjacent[0]
                    belongs_to._add_pixel(coord, intensity)
                else:
                    # Create a branch
                    belongs_to = Branch(adjacent, coord, intensity, idx=next_idx())
                    # Add branch to overall list
                    nodes[belongs_to.idx] = belongs_to
                
                # Set absolute index of pixel in index map
                self.index_map[coord] = belongs_to.idx

                # Add all insignificant leaves in 'merge' to the same object as this pixel:
                for m in merge:
                    # print "Merging leaf %i onto leaf %i" % (i, idx)
                    # Remove leaf
                    nodes.pop(m.idx)
                    # Merge the insignificant node that this pixel now belongs to:
                    belongs_to._merge(m)
                    # Update index map
                    m.add_footprint(self.index_map, belongs_to.idx)

        if verbose:
            progress_bar.progress = 100 # Done
            progress_bar.show_progress()
            print("") # newline

        # Create trunk from objects with no ancestors
        self.trunk = [node for node in nodes.itervalues() if node.parent is None]
        
        # Remove orphan leaves that aren't large enough        
        leaves_in_trunk = [node for node in self.trunk if type(node) == Leaf]
        for leaf in leaves_in_trunk:
            if (len(leaf.f) < min_npix or leaf.fmax - leaf.fmin < min_delta):
                # This leaf is an orphan, so remove all references to it:
                nodes.pop(leaf.idx)
                self.trunk.remove(leaf)
                leaf.add_footprint(self.index_map, 0)
        
        # To make the node.level property fast, we ensure all the nodes in the
        # trunk have their level cached as "0"
        for node in self.trunk:
            node._level = 0 # See the definition of level() in components.py
        
        # Save a list of all nodes accessible by ID
        self.nodes_dict = nodes 
        
        # Return the newly-created dendrogram:
        return self

    @staticmethod
    def _io_module_for(format):
        " Helper for load_from() and save_to() "
        try:
            # Load the module astrodendro.io.[format]
            return getattr(__import__('astrodendro.io.'+format).io, format)
        except ImportError:
            raise ValueError("Invalid/unknown format: '{}'".format(format))

    @staticmethod
    def load_from(filename, format="autodetect"):
        if format == "autodetect":
            format = filename.rsplit('.',1)[-1].lower()
        io_module = Dendrogram._io_module_for(format)
        return io_module.dendro_import(filename)
    
    def save_to(self, filename, format="autodetect"):
        if format == "autodetect":
            format = filename.rsplit('.',1)[-1].lower()
        io_module = Dendrogram._io_module_for(format)
        io_module.dendro_export(self, filename)
    
    @property
    def all_nodes(self):
        " Return a flattened iterable containing all nodes in the dendrogram "
        return self.nodes_dict.itervalues()
    
    def get_leaves(self):
        " Return a flattened list of all leaves in the dendrogram "
        return [i for i in self.nodes_dict.itervalues() if type(i) == Leaf]

    def to_newick(self):
        return "(%s);" % ','.join([node.newick for node in self.trunk])

    def node_at(self, coords):
        " Get the node at the given pixel coordinate, or None "
        idx = self.index_map[coords]
        if idx:
            return self.nodes_dict[idx]
        return None

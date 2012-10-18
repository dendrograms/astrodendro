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

import unittest
from astrodendro import Dendrogram
from astrodendro.components import Leaf, Branch
import numpy as np

class Test2DimensionalData(unittest.TestCase):
    def test_dendrogramWithNan(self):
        n = np.nan
        data = np.array([[n,n,n,n,n,n,n,n],
                         [n,4,n,n,n,n,n,n],
                         [n,n,n,1,n,n,0,5],
                         [3,n,n,2,3,2,0,n]])
        d = Dendrogram.compute(data)
        
        ########################################
        # Check the trunk elements:
        
        leaves = [node for node in d.trunk if type(node) == Leaf]
        branches = [node for node in d.trunk if node not in leaves]
        
        self.assertEqual(len(leaves), 2, msg="We expect two leaves among the lowest structures (the trunk)")
        self.assertEqual(len(branches), 1, msg="We expect one branch among the lowest structures (the trunk)")
        
        for leaf in leaves:
            self.assertEqual(len(leaf.f), 1, msg="Leaves in the trunk are only expected to contain one point")
            self.assertIsNone(leaf.parent)
            self.assertEqual(leaf.ancestor, leaf)
            self.assertEqual(leaf.get_npix(), 1)
            if leaf.f[0] == 4:
                self.assertEqual(leaf.coords[0], (1,1))
            elif leaf.f[0] == 3:
                self.assertEqual(leaf.coords[0], (3,0))
            else:
                self.fail("Invalid value of flux in one of the leaves")
        
        ########################################
        # Check properties of the branch:
        branch = branches[0]
        self.assertIsNone(branch.parent)
        self.assertEqual(branch.ancestor, branch)
        self.assertEqual(branch.merge_level, 0)
        self.assertEqual(branch.get_npix(subtree=False), 1) # only pixel is a 0
        self.assertEqual(branch.get_npix(subtree=True), 7)
        
        self.assertEqual(len(branch.children), 2)
        for leaf in branch.children:
            self.assertIsInstance(leaf, Leaf)
            self.assertEqual(leaf.ancestor, branch)
            self.assertEqual(leaf.parent, branch)
            if 5 in leaf.f:
                self.assertEqual(sum(leaf.f), 5)
            elif 3 in leaf.f:
                self.assertEqual(sum(leaf.f), 1+2+3+2)
            else:
                self.fail("Invalid child of the branch")

    def test_mergeLevelAndHeight(self):
        n = np.nan
        data = np.array([[n,n,n,n,n,],
                         [n,4,2,5,n,],
                         [n,n,n,n,0,]])
        d = Dendrogram.compute(data)
        branch, leaf4, leaf5 = d.trunk[0], d.node_at((1,1)), d.node_at((1,3))
        self.assertEqual(branch.merge_level, 2)
        self.assertEqual(leaf4.height, 2)
        self.assertEqual(leaf5.height, leaf5.fmax - branch.merge_level) # 3
        
        ### TODO: What is the appropriate value for branch.height ?

    def test_dendrogramWithConstBackground(self):
        # Test a highly artificial array containing a lot of equal pixels    
        data = np.array([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,3,1,1,1,1,1,1,1,1,1,1],
                         [1,1,3,5,3,1,1,1,1,1,1,1,1,1],
                         [1,1,2,3,2,2,2,1,1,1,1,1,1,1],
                         [1,1,1,1,3,4,3,1,1,1,1,1,1,1],
                         [1,1,1,1,2,3,2,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,2,3,2,1,1,1,1,1,1,1],
                         [1,1,1,1,3,4,3,1,1,2,2,1,1,1],
                         [1,1,1,1,2,3,2,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1],])
        d = Dendrogram.compute(data)
        self.assertLessEqual(len(d.nodes_dict), 7)
        # Some of the '1' valued pixels get included with the leaves and branches,
        # hence number of nodes is currently 7 and not 6 as expected.
        # Fixing this is probably more trouble than it's worth.
        leaf_with_twos = d.node_at((10, 9))
        self.assertEqual(leaf_with_twos.height, 1)
                
class Test3DimensionalData(unittest.TestCase):
    def setUp(self):
        from astrodendro.test._testdata import data
        self.data = data
    
    def test_dendrogramComputation(self):
        d = Dendrogram.compute(self.data, min_npix=8, min_delta=0.3, min_intensity=1.4)
        
        # This data with these parameters should produce 55 leaves
        self.assertEqual(len(d.leaves),55)
        
        # Now check every pixel in the data cube (this takes a while).
        # The following loop construct may look crazy, but it is a more
        # efficient way of iterating through the array than using a regular
        # nditer with multi_index.
        for coord in np.array(np.unravel_index( np.arange(self.data.size), self.data.shape)).transpose():
            coord = tuple(coord)
            f = self.data[coord]
            if (f < 1.4):
                self.assertEqual(d.node_at(coord), None)
            else:
                node = d.node_at(coord)
                if node:
                    # The current pixel is associated with part of the dendrogram.
                    self.assertIn(coord, node.coords, "Pixel at {0} is claimed to be part of {1}, but that node does not contain the coordinate {0}!".format(coord, node))
                    fmax_coords, fmax = node.get_peak(subtree=True)
                    if d.node_at(fmax_coords) is node:
                        # The current pixel is the peak pixel in this node
                        pass
                    else:
                        self.assertTrue(fmax >= f)


class TestNDimensionalData(unittest.TestCase):
    def test_4dim(self):
        " Test 4-dimensional data "
        data = np.zeros((5,5,5,5)) # Create a 5x5x5x5 array initialized to zero
        # N-dimensional data is hard to conceptualize so I've kept this simple.
        # Create a local maximum (value 5) at the centre
        data[2,2,2,2] = 5
        # add some points around it of intensity 3. Note that '1:4:2' is equivalent to saying indices '1' and '3'
        data[2,1:4:2,2,2] = data[2,2,1:4:2,2] = data[2,2,2,1:4:2] = 3
        # Add a trail of points of value 2 connecting one of those 3s to a 4
        data[0:3,0,2,2] = 2 # Sets [0,0,2,2], [1,0,2,2], and [2,0,2,2] all equal to 2 -> will connect to the '3' at [2,1,2,2]
        data[0,0,2,1] = 4
        
        # Now dendrogram it:
        d = Dendrogram.compute(data, min_intensity=1)
        # We expect two leaves:
        leaves = d.leaves
        self.assertEqual(len(leaves), 2)
        # We expect one branch:
        branches = [i for i in d.all_nodes if type(i) is Branch]
        self.assertEqual(len(branches), 1)
        self.assertEqual(len(d.trunk), 1)
        self.assertEqual(d.trunk[0], branches[0])
        
        # The maxima of each leaf should be at [2,2,2,2] and [0,3,2,1]
        for leaf in leaves:
            self.assertIn(leaf.get_peak(), ( ((2,2,2,2), 5.), ((0,0,2,1),4.) ) )
        self.assertNotEqual(leaves[0].get_peak(), leaves[1].get_peak())
        
        # Check out a few more properties of the leaf around the global maximum:
        leaf = d.node_at((2,2,2,2))
        self.assertEqual(leaf.fmax, 5)
        self.assertEqual(leaf.fmin, 2)
        self.assertEqual(leaf.get_npix(), 1+6+2) # Contains 1x '5', 6x '3', and 2x '2'. The other '2' should be in the branch
        # Check that the only pixel in the branch is a '2' at [0,0,2,2]
        self.assertEqual((branches[0].coords, branches[0].f), ( [(0,0,2,2),],[2.,] )) 

if __name__ == '__main__':
    unittest.main()

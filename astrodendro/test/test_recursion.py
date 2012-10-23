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

import sys

import numpy as np

from .. import Dendrogram


class TestRecursionLimit(object):
    """
    Test that we can efficiently compute deep dendrogram trees
    without hitting the recursion limit.
    Note: plot() uses recursion but we should be able to *compute*
    dendrograms without using deep recursion, even if we aren't
    yet able to plot them without using recursion.
    """
    def setup_method(self, method):
        self._oldlimit = sys.getrecursionlimit()
        sys.setrecursionlimit(100)  # Reduce recursion limit dramatically (default is 1000)
        size = 10000  # number of leaves desired in the dendrogram
        data1 = np.arange(size * 2)  # first row
        data2 = np.arange(size * 2)  # second row
        data2[::2] += 2
        data1[-1] = 0  # set the last pixel in the first row to zero, to trigger a deep ancestor search
        self.data = np.vstack((data1, data2))
        self.size = size
        # self.data now looks like this:
        # [[ 0, 1, 2, 3, 4, 5, ...],
        #  [ 2, 1, 4, 3, 6, 5, ...]]
        # Notice every second column has a local maximum
        # so there are [size] local maxima in the array

    def test_compute(self):
        d = Dendrogram.compute(self.data)
        assert len(d.leaves) == self.size, "We expect {n} leaves, not {a}.".format(n=self.size, a=len(d.leaves))

    def test_computing_level(self):
        d = Dendrogram.compute(self.data)

        # Now pick a node near the middle of the dendrogram:
        mid_node = d.node_at((0, self.size // 2))

        # Compute its level:
        sys.setrecursionlimit(100000)
        _ = mid_node.level

        # Now check the .level property of all nodes, in random order:
        import random
        nodes = random.sample(list(d.all_nodes), len(d.nodes_dict))
        for node in nodes:
            obj = node
            level = node.level
            while level > 0:
                obj = obj.parent
                level -= 1
            assert obj.parent == None
            assert obj.level == 0

    def teardown_method(self, method):
        sys.setrecursionlimit(self._oldlimit)

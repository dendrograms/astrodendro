# Licensed under an MIT open source license - see LICENSE

import sys

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')

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
        self._make_data(size)

    def _make_data(self, size):
        data1 = np.arange(size * 2)  # first row
        data2 = np.arange(size * 2)  # second row
        data2[::2] += 2
        data1[-1] = 0  # set the last pixel in the first row to zero, to trigger a deep ancestor search
        data = np.vstack((data1, data2))
        self.data = data
        self.size = size

        # result looks like this:
        # [[ 0, 1, 2, 3, 4, 5, ...],
        #  [ 2, 1, 4, 3, 6, 5, ...]]
        # Notice every second column has a local maximum
        # so there are [size] local maxima in the array

    def test_compute(self):
        d = Dendrogram.compute(self.data)
        assert len(d.leaves) == self.size, "We expect {n} leaves, not {a}.".format(n=self.size, a=len(d.leaves))

    def test_computing_level(self):
        d = Dendrogram.compute(self.data)

        # Now pick a structure near the middle of the dendrogram:
        mid_structure = d.structure_at((0, self.size // 2))

        # Compute its level:
        sys.setrecursionlimit(100000)
        _ = mid_structure.level

        # Check that .level satisfies the recurrence relation:
        # 0 if root else parent.level + 1
        for structure in d.all_structures:
            if structure.parent is None:
                assert structure.level == 0
            else:
                assert structure.level == structure.parent.level + 1


    def test_plot(self):
        sys.setrecursionlimit(self._oldlimit)
        ax = plt.gca()
        sys.setrecursionlimit(150)

        d = Dendrogram.compute(self.data)
        p = d.plotter()
        p.plot_tree(ax)

    def teardown_method(self, method):
        sys.setrecursionlimit(self._oldlimit)

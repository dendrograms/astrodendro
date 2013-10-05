# Licensed under an MIT open source license - see LICENSE

import numpy as np
from ..dendrogram import Dendrogram, TreeIndex

def assert_permuted_fancyindex(x, y):
    """ Assert that two fancy indices (tuples of integer ndarrays)
    are permutations of each other
    """
    if not isinstance(x, tuple) or not(isinstance(x[0], np.ndarray)):
        raise TypeError("First argument not a fancy index: %s" % x)

    if not isinstance(y, tuple) or not(isinstance(y[0], np.ndarray)):
        raise TypeError("Second argument not a fancy index: %s" % y)

    dtype = [('%i' % i, 'i') for i in range(len(x))]
    x = np.array(list(zip(*x)), dtype=dtype)
    y = np.array(list(zip(*y)), dtype=dtype)
    np.testing.assert_array_equal(np.sort(x),
                                  np.sort(y))

def assert_identical_fancyindex(x, y):
    for xx, yy in zip(x, y):
        np.testing.assert_array_equal(xx, yy)


class TestIndex(object):

    def setup_method(self, method):
        pass

    def assert_valid_index(self, d, index):
        """Assert that a dendrogram index is correct"""

        #subtree=False is a permutation of np.where(index_map == x)
        for s in d.all_structures:
            ind = index.indices(s.idx, subtree=False)
            expected = np.where(d.index_map == s.idx)
            assert_permuted_fancyindex(ind, expected)

        #subtree=True is the same, but includes descendents
        for s in d.all_structures:
            ind = index.indices(s.idx, subtree=True)
            expected = [np.where(d.index_map == ss.idx) for
                        ss in [s] + s.descendants]
            expected = tuple(np.hstack(i) for i in zip(*expected))
            assert_permuted_fancyindex(ind, expected)

    def test_single_trunk(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 3, 5],
                         [1, 1, 1, 4]])
        d = Dendrogram.compute(data)
        assert len(d) == 3
        self.assert_valid_index(d, TreeIndex(d))

    def test_two_trunk(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 1, 5],
                         [1, 1, 1, 4]])
        d = Dendrogram.compute(data, min_value=2)
        assert len(d.trunk) == 2
        self.assert_valid_index(d, TreeIndex(d))

    def test_single_structure(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 1, 1],
                         [1, 1, 1, 1]])
        d = Dendrogram.compute(data)
        assert len(d) == 1
        self.assert_valid_index(d, TreeIndex(d))

    def test_cube(self):
        np.random.seed(42)
        data = np.random.random((5, 5, 5))

        d = Dendrogram.compute(data)
        self.assert_valid_index(d, TreeIndex(d))

    def test_1d(self):
        np.random.seed(42)
        data = np.random.random(5)

        d = Dendrogram.compute(data)
        assert len(d) > 0
        self.assert_valid_index(d, TreeIndex(d))

    def test_4d(self):
        np.random.seed(42)
        data = np.random.random((3, 3, 3, 3))

        d = Dendrogram.compute(data)
        assert len(d) > 0
        self.assert_valid_index(d, TreeIndex(d))

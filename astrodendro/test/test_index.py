import numpy as np
from ..dendrogram import Dendrogram, TreeIndex


class TestIndex(object):

    def setup_method(self, method):
        pass

    def assert_equivalent_indices(self, x, y):
        """Assert that two multidimensional indices
        are permutations of each other"""
        dtype = [('%i' % i, 'i') for i in range(len(x))]
        x = np.array(zip(*x), dtype=dtype)
        y = np.array(zip(*y), dtype=dtype)
        np.testing.assert_array_equal(np.sort(x),
                                      np.sort(y))

    def assert_valid_index(self, d, index):
        """Assert that a dendrogram index is correct"""

        #subtree=False is a permutation of np.where(index_map == x)
        for s in d.all_nodes:
            ind = index.indices(s.idx)
            expected = np.where(d.index_map == s.idx)
            self.assert_equivalent_indices(ind, expected)

        #subtree=True is the same, but includes descendents
        for s in d.all_nodes:
            ind = index.indices(s.idx, subtree=True)
            expected = [np.where(d.index_map == ss.idx) for
                        ss in [s] + s.descendants]
            expected = tuple(np.hstack(i) for i in zip(*expected))
            self.assert_equivalent_indices(ind, expected)

    def test_single_trunk(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 3, 5],
                         [1, 1, 1, 4]])
        d = Dendrogram.compute(data)
        assert len(d.nodes_dict) == 3
        self.assert_valid_index(d, TreeIndex(d))

    def test_two_trunk(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 1, 5],
                         [1, 1, 1, 4]])
        d = Dendrogram.compute(data, min_data_value=2)
        assert len(d.trunk) == 2
        self.assert_valid_index(d, TreeIndex(d))

    def test_single_structure(self):
        data = np.array([[1, 1, 1, 1],
                         [1, 5, 1, 1],
                         [1, 4, 1, 1],
                         [1, 1, 1, 1]])
        d = Dendrogram.compute(data)
        assert len(d.nodes_dict) == 1
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
        assert len(d.nodes_dict) > 0
        self.assert_valid_index(d, TreeIndex(d))

    def test_4d(self):
        np.random.seed(42)
        data = np.random.random((3, 3, 3, 3))

        d = Dendrogram.compute(data)
        assert len(d.nodes_dict) > 0
        self.assert_valid_index(d, TreeIndex(d))

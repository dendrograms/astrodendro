import numpy as np

from ..pruning import *
from ..structure import Structure
from ..dendrogram import Dendrogram
from .test_index import assert_permuted_fancyindex


data = np.array([[0, 0, 0, 0],
                 [0, 2, 1, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 0]])
st = Structure(zip(*np.where(data)), data[np.where(data)])

def test_npix():
    assert min_npix(3)(st)
    assert not min_npix(4)(st)

def test_min_delta():
    assert min_delta(1)(st)
    assert not min_delta(1.1)(st)
    assert min_delta(1.1)(st, value=0)


def test_min_peak():
    assert min_peak(2)(st)
    assert not min_peak(2.1)(st)

def test_min_sum():
    assert min_sum(4)(st)
    assert not min_sum(4.1)(st)

def test_contains_seeds():
    assert contains_seeds(np.where(data == 2))(st)
    assert contains_seeds(np.where(data == 1))(st)
    assert contains_seeds(([1], [1]))(st)
    assert not contains_seeds(np.where(data == 0))(st)

def test_all_true():
    c1 = min_npix(3)
    c2 = min_peak(2)
    c3 = min_npix(4)
    assert all_true((c1, c2))(st)
    assert not all_true((c1, c2, c3))(st)

def test_multi_ravel():
    from ..pruning import _ravel_multi_index

    x = _ravel_multi_index([[0, 1], [0, 1]], [3, 3])
    np.testing.assert_array_equal(x, [0, 4])

    x = _ravel_multi_index([[1, 0], [0, 1]], [4, 2])
    np.testing.assert_array_equal(x, [2, 1])

    x = _ravel_multi_index([[0], [9]], [5, 3], mode='clip')
    np.testing.assert_array_equal(x, [2])

    x = _ravel_multi_index([[0], [9]], [5, 3], mode='wrap')
    np.testing.assert_array_equal(x, [0])

def compare_dendrograms(d1, d2):

    assert len(d1) == len(d2)
    assert (np.sort([leaf.vmax for leaf in d1.all_structures]) == np.sort([leaf.vmax for leaf in d2.all_structures])).all()
    assert (np.sort([leaf.vmin for leaf in d1.all_structures]) == np.sort([leaf.vmin for leaf in d2.all_structures])).all()

    for s in d1.all_structures:
         ind1 = np.where(d1.index_map == s.idx)
         idx2 = d2.index_map[ind1][0]
         ind2 = np.where(d2.index_map == idx2)
         assert_permuted_fancyindex(ind1, ind2)

class TestPostPruning(object):

    def setup_class(self):
        from ._testdata import data
        self.data = data

    def test_min_delta(self):
        d1 = Dendrogram.compute(self.data, min_delta=1.0, min_npix=8, min_value=1.4)

        d2 = Dendrogram.compute(self.data, min_npix=8, min_value=1.4)
        d2.prune(min_delta=1.0)

        compare_dendrograms(d1, d2)

    def test_min_npix(self):
        d1 = Dendrogram.compute(self.data, min_npix=8, min_delta=0.3, min_value=1.4)

        d2 = Dendrogram.compute(self.data, min_delta=0.3, min_value=1.4)
        d2.prune(min_npix=8)

        compare_dendrograms(d1, d2)

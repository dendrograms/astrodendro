# The tests here ensure that the Structure class behaves as expected and has
# the correct interface. This is done by setting up a few simple examples and
# ensuring that all the properties and methods are behaving as expected.

from itertools import repeat

import pytest
import numpy as np
from numpy.testing import assert_allclose
from ..structure import Structure


@pytest.mark.parametrize('index', [(0,), (1, 3), (4, 5, 9)])
def test_init_leaf_scalar(index):

    s = Structure(index, 1.5)

    # Properties
    assert s.idx is None
    assert s.is_leaf
    assert not s.is_branch
    assert np.all(s.indices == np.array([index]))
    assert np.all(s.indices_all == s.indices)
    assert np.all(s.values == np.array([1.5]))
    assert np.all(s.values_all == s.values)
    assert s.vmin == 1.5
    assert s.vmax == 1.5
    assert s.height == 0.
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == []
    assert s.descendants == []

    # Methods
    for subtree in [False, True]:
        assert s.get_npix(subtree=subtree) == 1
        assert s.get_peak(subtree=subtree) == (index, 1.5)

    # Footprint
    array = np.zeros([20 for i in range(len(index))])
    s.fill_footprint(array, level=2)
    assert array[index] == 2
    assert np.sum(array) == 2.


@pytest.mark.parametrize('index', [[(0,), (1,), (2,)],
                                   [(1, 3), (2, 2), (4, 1)],
                                   [(4, 5, 9), (3, 2, 1), (6, 7, 8)]])
def test_init_leaf_list(index):

    s = Structure(index, [3.1, 4.2, 5.3])

    # Properties
    assert s.idx is None
    assert s.is_leaf
    assert not s.is_branch
    assert np.all(s.indices == np.array(index))
    assert np.all(s.indices_all == s.indices)
    assert np.all(s.values == np.array([3.1, 4.2, 5.3]))
    assert np.all(s.values_all == s.values)
    assert s.vmin == 3.1
    assert s.vmax == 5.3
    assert_allclose(s.height, 2.2)
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == []
    assert s.descendants == []

    # Methods
    for subtree in [False, True]:
        assert s.get_npix(subtree=subtree) == 3
        assert s.get_peak(subtree=subtree) == (index[2], 5.3)

    # Footprint
    array = np.zeros([20 for i in range(len(index[0]))])
    s.fill_footprint(array, level=2)
    for i in index:
        print(i, array[i])
        assert array[i] == 2
    assert np.sum(array) == 6.


@pytest.mark.parametrize('index', [(0,), (1, 3), (4, 5, 9)])
def test_init_branch_scalar(index):

    leaf_index = tuple([10 for i in range(len(index))])
    leaf = Structure(leaf_index, 20.)

    s = Structure(index, 1.5, children=[leaf])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch
    assert np.all(s.indices == [index])
    assert np.all(s.indices_all == np.array([index] + leaf.indices))
    assert np.all(s.values == np.array([1.5]))
    assert np.all(s.values_all == np.array([1.5] + leaf.values))
    assert s.vmin == 1.5
    assert s.vmax == 1.5
    assert s.height == 0.
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == [leaf]
    assert s.descendants == [leaf]

    # Leaf properties
    assert leaf.level == 1
    assert leaf.ancestor is s
    assert leaf.parent is s
    assert leaf.children == []
    assert leaf.descendants == []

    # Methods
    assert s.get_npix(subtree=False) == 1
    assert s.get_peak(subtree=False) == (index, 1.5)
    assert s.get_npix(subtree=True) == 2
    assert s.get_peak(subtree=True) == (leaf_index, 20.)

    # Footprint
    array = np.zeros([20 for i in range(len(index))])
    s.fill_footprint(array, level=2)
    assert array[index] == 2.
    assert array[leaf_index] == 3.
    assert np.sum(array) == 5.


@pytest.mark.parametrize('index', [[(0,), (1,), (2,)],
                                   [(1, 3), (2, 2), (4, 1)],
                                   [(4, 5, 9), (3, 2, 1), (6, 7, 8)]])
def test_init_branch_list(index):

    ndim = len(index[0])

    leaf_index = tuple([10 for i in range(ndim)])

    leaf = Structure(leaf_index, 20.)

    s = Structure(index, [3.1, 4.2, 5.3], children=[leaf])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch
    assert np.all(s.indices == np.array(index))
    print index, leaf.indices
    print index + leaf.indices
    assert np.all(s.indices_all == np.array(index + leaf.indices))
    assert np.all(s.values == np.array([3.1, 4.2, 5.3]))
    assert np.all(s.values_all == np.array([3.1, 4.2, 5.3] + leaf.values))
    assert s.vmin == 3.1
    assert s.vmax == 5.3
    assert_allclose(s.height, 2.2)
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == [leaf]
    assert s.descendants == [leaf]

    # Leaf properties
    assert leaf.level == 1
    assert leaf.ancestor is s
    assert leaf.parent is s
    assert leaf.children == []
    assert leaf.descendants == []

    # Methods
    assert s.get_npix(subtree=False) == 3
    assert s.get_peak(subtree=False) == (index[2], 5.3)
    assert s.get_npix(subtree=True) == 4
    assert s.get_peak(subtree=True) == (leaf_index, 20.)

    # Footprint
    array = np.zeros([20 for i in range(ndim)])
    s.fill_footprint(array, level=2)
    for i in index:
        assert array[i] == 2.
    assert array[leaf_index] == 3.
    assert np.sum(array) == 9.

# TODO: add newick tests

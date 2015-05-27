# Licensed under an MIT open source license - see LICENSE

# The tests here ensure that the Structure class behaves as expected and has
# the correct interface. This is done by setting up a few simple examples and
# ensuring that all the properties and methods are behaving as expected.

import pytest
import numpy as np
from numpy.testing import assert_allclose
from ..structure import Structure
from .test_index import assert_identical_fancyindex


@pytest.mark.parametrize('index', [(0,), (1, 3), (4, 5, 9)])
def test_init_leaf_scalar(index):

    s = Structure(index, 1.5)

    # Properties
    assert s.idx is None
    assert s.is_leaf
    assert not s.is_branch
    assert_identical_fancyindex(s.indices(subtree=False),
                               tuple(np.atleast_1d(i) for i in index))
    assert np.all(s.indices() == s.indices(subtree=True))
    assert np.all(s.values(subtree=False) == np.array([1.5]))
    assert np.all(s.values(subtree=True) == s.values(subtree=False))
    assert s.vmin == 1.5
    assert s.vmax == 1.5
    assert s.height == 1.5
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
    s._fill_footprint(array, level=2)
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
    indices = tuple(np.atleast_1d(i) for i in zip(*index))
    assert_identical_fancyindex(s.indices(subtree=False), indices)
    assert_identical_fancyindex(s.indices(subtree=True), indices)
    assert np.all(s.values(subtree=False) == np.array([3.1, 4.2, 5.3]))
    assert np.all(s.values(subtree=True) == s.values(subtree=False))
    assert s.vmin == 3.1
    assert s.vmax == 5.3
    assert s.height == 5.3
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
    s._fill_footprint(array, level=2)
    for i in index:
        print(i, array[i])
        assert array[i] == 2
    assert np.sum(array) == 6.


@pytest.mark.parametrize('index', [(0,), (1, 3), (4, 5, 9)])
def test_init_branch_scalar(index):

    leaf_index = tuple([10 for i in range(len(index))])
    leaf = Structure(leaf_index, 20.)
    leaf_indices = leaf.indices(subtree=False)

    s = Structure(index, 1.5, children=[leaf])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch
    indices = tuple(np.atleast_1d(i) for i in index)
    indices_all = tuple(np.hstack(a)
                        for a in zip(indices, leaf_indices))
    assert_identical_fancyindex(s.indices(subtree=False), indices)
    assert_identical_fancyindex(s.indices(subtree=True), indices_all)
    assert np.all(s.values(subtree=False) == np.array([1.5]))
    assert np.all(s.values(subtree=True) == np.hstack(([1.5], leaf.values(subtree=False))))
    assert s.vmin == 1.5
    assert s.vmax == 1.5
    assert s.height == 20.
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
    s._fill_footprint(array, level=2)
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
    leaf_indices = leaf.indices(subtree=False)

    s = Structure(index, [3.1, 4.2, 5.3], children=[leaf])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch
    indices = tuple(np.atleast_1d(i) for i in zip(*index))
    indices_all = tuple(np.hstack(a) for a in
                        zip(indices, leaf_indices))

    assert_identical_fancyindex(s.indices(subtree=False), indices)
    assert_identical_fancyindex(s.indices(subtree=True), indices_all)

    assert np.all(s.values(subtree=False) == np.array([3.1, 4.2, 5.3]))
    assert np.all(s.values(subtree=True) == np.hstack(([3.1, 4.2, 5.3], leaf.values(subtree=False))))
    assert s.vmin == 3.1
    assert s.vmax == 5.3
    assert s.height == 20.
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
    s._fill_footprint(array, level=2)
    for i in index:
        assert array[i] == 2.
    assert array[leaf_index] == 3.
    assert np.sum(array) == 9.


@pytest.mark.parametrize('index', [(0,), (1, 3), (4, 5, 9)])
def test_init_branch_scalar_3_level(index):

    leaf_index = tuple([10 for i in range(len(index))])
    leaf = Structure(leaf_index, 20.)
    leaf_indices = leaf.indices(subtree=False)

    branch_index = tuple([9 for i in range(len(index))])
    branch = Structure(branch_index, 15., children=[leaf])
    branch_indices = branch.indices(subtree=False)

    s = Structure(index, 1.5, children=[branch])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch

    indices = tuple(np.atleast_1d(i) for i in index)
    indices_all = tuple(np.hstack(a)
                        for a in zip(indices, branch_indices, leaf_indices))
    assert_identical_fancyindex(s.indices(subtree=False), indices)
    assert_identical_fancyindex(s.indices(subtree=True), indices_all)

    assert np.all(s.values(subtree=False) == np.array([1.5]))
    assert np.all(s.values(subtree=True) == np.hstack((s.values(subtree=False), branch.values(subtree=False), leaf.values(subtree=False))))
    assert s.vmin == 1.5
    assert s.vmax == 1.5
    assert s.height == 15.
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == [branch]
    assert set(s.descendants) == set([branch, leaf])  # order should not matter

    # Branch properties
    assert branch.level == 1
    assert branch.ancestor is s
    assert branch.parent is s
    assert branch.children == [leaf]
    assert branch.descendants == [leaf]

    # Leaf properties
    assert leaf.level == 2
    assert leaf.ancestor is s
    assert leaf.parent is branch
    assert leaf.children == []
    assert leaf.descendants == []

    # Methods
    assert s.get_npix(subtree=False) == 1
    assert s.get_peak(subtree=False) == (index, 1.5)
    assert s.get_npix(subtree=True) == 3
    assert s.get_peak(subtree=True) == (leaf_index, 20.)

    # Footprint
    array = np.zeros([20 for i in range(len(index))])
    s._fill_footprint(array, level=2)
    assert array[index] == 2.
    assert array[branch_index] == 3.
    assert array[leaf_index] == 4.
    assert np.sum(array) == 9.


@pytest.mark.parametrize('index', [[(0,), (1,), (2,)],
                                   [(1, 3), (2, 2), (4, 1)],
                                   [(4, 5, 9), (3, 2, 1), (6, 7, 8)]])
def test_init_branch_list_3_level(index):

    ndim = len(index[0])

    leaf_index = tuple([10 for i in range(ndim)])
    leaf = Structure(leaf_index, 20.)
    leaf_indices = leaf.indices(subtree=False)

    branch_index = tuple([9 for i in range(ndim)])
    branch = Structure(branch_index, 15., children=[leaf])
    branch_indices = branch.indices(subtree=False)

    s = Structure(index, [3.1, 4.2, 5.3], children=[branch])

    # Properties
    assert s.idx is None
    assert not s.is_leaf
    assert s.is_branch

    indices = tuple(np.atleast_1d(i) for i in zip(*index))
    indices_all = tuple(np.hstack(a) for a in
                        zip(indices, branch_indices, leaf_indices))
    assert_identical_fancyindex(s.indices(subtree=False), indices)
    assert_identical_fancyindex(s.indices(subtree=True), indices_all)

    assert np.all(s.values(subtree=False) == np.array([3.1, 4.2, 5.3]))
    assert np.all(s.values(subtree=True) == np.hstack((s.values(subtree=False), branch.values(subtree=False), leaf.values(subtree=False))))
    assert s.vmin == 3.1
    assert s.vmax == 5.3
    assert s.height == 15.
    assert s.level == 0
    assert s.ancestor is s
    assert s.parent is None
    assert s.children == [branch]
    assert s.descendants == [branch, leaf]

    # Branch properties
    assert branch.level == 1
    assert branch.ancestor is s
    assert branch.parent is s
    assert branch.children == [leaf]
    assert branch.descendants == [leaf]

    # Leaf properties
    assert leaf.level == 2
    assert leaf.ancestor is s
    assert leaf.parent is branch
    assert leaf.children == []
    assert leaf.descendants == []

    # Methods
    assert s.get_npix(subtree=False) == 3
    assert s.get_peak(subtree=False) == (index[2], 5.3)
    assert s.get_npix(subtree=True) == 5
    assert s.get_peak(subtree=True) == (leaf_index, 20.)

    # Footprint
    array = np.zeros([20 for i in range(ndim)])
    s._fill_footprint(array, level=2)
    for i in index:
        assert array[i] == 2.
    assert array[branch_index] == 3.
    assert array[leaf_index] == 4.
    assert np.sum(array) == 13.


def test_add_pixel():

    s = Structure(1, 10.)

    assert s.get_npix() == 1
    assert s.get_peak() == (1, 10)
    assert s.vmin == 10.
    assert s.vmax == 10.

    s._add_pixel(2, 8.)

    assert s.get_npix() == 2
    assert s.get_peak() == (1, 10)
    assert s.vmin == 8.
    assert s.vmax == 10.

    s._add_pixel(3, 12.)

    assert s.get_npix() == 3
    assert s.get_peak() == (3, 12.)
    assert s.vmin == 8.
    assert s.vmax == 12.

def test_sorted_leaves():
    l1 = Structure(1, 10., idx=1)
    l2 = Structure(2, 8., idx=2)
    s = Structure(3, 5., children=[l1, l2], idx=3)
    assert s.sorted_leaves() == [l2, l1]
    assert s.sorted_leaves(reverse=True) == [l1, l2]
    def key(x):
        return x.idx
    assert s.sorted_leaves(sort_key=key) == [l1, l2]

    s2 = Structure(4, 3., children=[s], idx=4)
    assert s2.sorted_leaves() == [l2, l1]
    assert s2.sorted_leaves(subtree=False) == []

# TODO: add newick tests

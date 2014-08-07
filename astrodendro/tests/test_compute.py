# Licensed under an MIT open source license - see LICENSE

import numpy as np
import pytest

from .. import Dendrogram, periodic_neighbours, Structure


class Test2DimensionalData(object):

    def test_dendrogramWithNan(self):

        n = np.nan
        data = np.array([[n, n, n, n, n, n, n, n],
                         [n, 4, n, n, n, n, n, n],
                         [n, n, n, 1, n, n, 0, 5],
                         [3, n, n, 2, 3, 2, 0, n]])
        d = Dendrogram.compute(data)

        ########################################
        # Check the trunk elements:

        leaves = [structure for structure in d.trunk if structure.is_leaf]
        branches = [structure for structure in d.trunk if structure not in leaves]

        assert len(leaves) == 2, "We expect two leaves among the lowest structures (the trunk)"
        assert len(branches) == 1, "We expect one branch among the lowest structures (the trunk)"

        for leaf in leaves:
            assert len(leaf.values(subtree=False)) == 1, "Leaves in the trunk are only expected to contain one point"
            assert leaf.parent is None
            assert leaf.ancestor == leaf
            assert leaf.get_npix() == 1
            if leaf.values(subtree=False)[0] == 4:
                assert list(zip(*leaf.indices(subtree=False)))[0] == (1, 1)
            elif leaf.values(subtree=False)[0] == 3:
                assert list(zip(*leaf.indices(subtree=False)))[0] == (3, 0)
            else:
                self.fail("Invalid value of flux in one of the leaves")

        ########################################
        # Check properties of the branch:
        branch = branches[0]
        assert branch.parent is None
        assert branch.ancestor == branch
        assert branch.get_npix(subtree=False) == 1  # only pixel is a 0
        assert branch.get_npix(subtree=True) == 7

        assert len(branch.children) == 2
        for leaf in branch.children:
            assert leaf.is_leaf
            assert leaf.ancestor == branch
            assert leaf.parent == branch
            if 5 in leaf.values(subtree=False):
                assert sum(leaf.values(subtree=False)) == 5
            elif 3 in leaf.values(subtree=False):
                assert sum(leaf.values(subtree=False)) == 1 + 2 + 3 + 2
            else:
                self.fail("Invalid child of the branch")

    def test_mergeLevelAndHeight(self):
        n = np.nan
        data = np.array([[n, n, n, n, n, ],
                         [n, 4, 2, 5, n, ],
                         [n, n, n, n, 0, ]])
        d = Dendrogram.compute(data)
        branch, leaf4, leaf5 = d.trunk[0], d.structure_at((1, 1)), d.structure_at((1, 3))
        assert leaf4.height == 4.
        assert leaf5.height == 5.
        assert branch.height == 4.

    def test_dendrogramWithConstBackground(self):
        # Test a highly artificial array containing a lot of equal pixels
        data = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 3, 5, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 3, 4, 3, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 3, 4, 3, 1, 1, 2, 2, 1, 1, 1],
                         [1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ])
        d = Dendrogram.compute(data)
        assert len(d) <= 7
        # Some of the '1' valued pixels get included with the leaves and branches,
        # hence number of structures is currently 7 and not 6 as expected.
        # Fixing this is probably more trouble than it's worth.
        leaf_with_twos = d.structure_at((10, 9))
        assert leaf_with_twos.height == 2

        # Check that all structures contain a reference to the dendrogram
        for structure in d:
            assert structure._dendrogram is d


class Test3DimensionalData(object):
    def setup_method(self, method):
        from ._testdata import data
        self.data = data

    def test_dendrogramComputation(self):
        d = Dendrogram.compute(self.data, min_npix=8, min_delta=0.3, min_value=1.4)

        # This data with these parameters should produce 55 leaves
        assert len(d.leaves) == 55

        # Now check every pixel in the data cube (this takes a while).
        st_map = -np.ones(self.data.shape, dtype=np.int)
        for st in d.all_structures:
            st_map[st.indices(subtree=False)] = st.idx

        #check that vmin/vmax/peak are correct
        for st in d.all_structures:
            assert st.vmin == self.data[st.indices(subtree=False)].min()
            assert st.vmax == self.data[st.indices(subtree=False)].max()
            pk_exp = self.data[st.indices(subtree=True)].max()

            ind, pk = st.get_peak(subtree=True)
            assert self.data[ind] == pk
            assert pk_exp == pk

        # The "right" way to do this is loop through indices,
        # and repeatedly call structure_at(). However, this is quite slow
        # structure_at is a thin wrapper around index_map,
        # and we compare index_map to st_map instead
        np.testing.assert_array_equal(st_map, d.index_map)

        # here, we test a few values of structure_at
        for coord in np.indices(self.data.shape).reshape(self.data.ndim, np.prod(self.data.shape)).transpose()[::100]:
            coord = tuple(coord)
            f = self.data[coord]
            structure = d.structure_at(coord)
            if structure is not None:
                assert structure.idx == st_map[coord], "Pixel at {0} is claimed to be part of {1}, but that structure does not contain the coordinate {0}!".format(coord, structure)
            else:
                assert st_map[coord] == -1


class TestNDimensionalData(object):
    def test_4dim(self):
        " Test 4-dimensional data "
        data = np.zeros((5, 5, 5, 5))  # Create a 5x5x5x5 array initialized to zero
        # N-dimensional data is hard to conceptualize so I've kept this simple.
        # Create a local maximum (value 5) at the centre
        data[2, 2, 2, 2] = 5
        # add some points around it with value 3. Note that '1:4:2' is equivalent to saying indices '1' and '3'
        data[2, 1:4:2, 2, 2] = data[2, 2, 1:4:2, 2] = data[2, 2, 2, 1:4:2] = 3
        # Add a trail of points of value 2 connecting one of those 3s to a 4
        data[0:3, 0, 2, 2] = 2  # Sets [0, 0, 2, 2], [1, 0, 2, 2], and [2, 0, 2, 2] all equal to 2 -> will connect to the '3' at [2, 1, 2, 2]
        data[0, 0, 2, 1] = 4

        # Now dendrogram it:
        d = Dendrogram.compute(data, min_value=1)
        # We expect two leaves:
        leaves = d.leaves
        assert len(leaves) == 2
        # We expect one branch:
        branches = [i for i in d.all_structures if i.is_branch]
        assert len(branches) == 1
        assert len(d.trunk) == 1
        assert d.trunk[0] == branches[0]

        # The maxima of each leaf should be at [2,2,2,2] and [0,3,2,1]
        for leaf in leaves:
            assert leaf.get_peak() in (((2, 2, 2, 2), 5.), ((0, 0, 2, 1), 4.))
        assert leaves[0].get_peak() != leaves[1].get_peak()

        # Check out a few more properties of the leaf around the global maximum:
        leaf = d.structure_at((2, 2, 2, 2))
        assert leaf.vmax == 5
        assert leaf.vmin == 2
        assert leaf.get_npix() == 1 + 6 + 2  # Contains 1x '5', 6x '3', and 2x '2'. The other '2' should be in the branch
        # Check that the only pixel in the branch is a '2' at [0,0,2,2]
        assert (list(zip(*branches[0].indices(subtree=False))), branches[0].values(subtree=False)) == ([(0, 0, 2, 2), ], [2., ])


def test_periodic():
    x = np.array([[0, 0, 0, 0, 0, ],
                 [1, 1, 0, 1, 1],
                 [0, 0, 0, 0, 0]])

    d = Dendrogram.compute(x, min_value=0.5,
                           neighbours=periodic_neighbours(1))
    expected = np.array([[-1, -1, -1, -1, -1],
                        [0, 0, -1, 0, 0],
                        [-1, -1, -1, -1, -1]])
    np.testing.assert_array_equal(d.index_map, expected)

def test_periodic_left():
    x = np.array([[1, 0, 0, 0, 0],
                  [1, 0, 0, 0, 1],
                  [1, 0, 0, 0, 0]])
    d = Dendrogram.compute(x, min_value=0.5,
                           neighbours=periodic_neighbours(1))
    expected = np.array([[0, -1, -1, -1, -1],
                         [0, -1, -1, -1, 0],
                         [0, -1, -1, -1, -1]])
    np.testing.assert_array_equal(d.index_map, expected)

def test_periodic_left_narrow():
    x = np.array([[0, 0, 0, 0, 0],
                  [1, 1, 0, 0, 1],
                  [0, 0, 0, 0, 0]])
    d = Dendrogram.compute(x, min_value=0.5,
                           neighbours=periodic_neighbours(1))
    expected = np.array([[-1, -1, -1, -1, -1],
                         [0, 0, -1, -1, 0],
                         [-1, -1, -1, -1, -1]])
    np.testing.assert_array_equal(d.index_map, expected)

def test_periodic_right():
    x = np.array([[0, 0, 0, 0, 1],
                  [1, 0, 0, 0, 1],
                  [0, 0, 0, 0, 1]])
    d = Dendrogram.compute(x, min_value=0.5,
                           neighbours=periodic_neighbours(1))
    expected = np.array([[-1, -1, -1, -1, 0],
                         [0, -1, -1, -1, 0],
                         [-1, -1, -1, -1, 0]])
    np.testing.assert_array_equal(d.index_map, expected)

def test_periodic_right_narrow():
    x = np.array([[0, 0, 0, 0, 0],
                  [1, 0, 0, 1, 1],
                  [0, 0, 0, 0, 0]])
    d = Dendrogram.compute(x, min_value=0.5,
                           neighbours=periodic_neighbours(1))
    expected = np.array([[-1, -1, -1, -1, -1],
                         [0, -1, -1, 0, 0],
                         [-1, -1, -1, -1, -1]])
    np.testing.assert_array_equal(d.index_map, expected)



from .build_benchmark import BENCHMARKS


@pytest.mark.parametrize(('filename'), BENCHMARKS.keys())
def test_benchmark(filename):
    from astropy.io import fits
    import os

    path = os.path.join(os.path.dirname(__file__),
                        'benchmark_data', filename)
    p = BENCHMARKS[filename]
    data = fits.getdata(path, 1)

    d1 = Dendrogram.compute(data, **p)
    d2 = Dendrogram.load_from(path)

    assert d1 == d2

    # Check that all structures contain a reference to the dendrogram
    for structure in d1:
        assert structure._dendrogram is d1


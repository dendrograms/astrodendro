# These tests ensure that the ``is_independent`` function is working correctly

import pytest
import numpy as np

from ..dendrogram import Dendrogram


class TestCustomMerge(object):

    def setup_class(self):
        self.data = np.array([0, 3.3, 5.5, 2.1, 1.0, 6.0, 4.4, 1.5, 4.1, 0.5])

    def test_reference(self):

        d = Dendrogram.compute(self.data)

        branches = [s for s in d.all_structures if s.is_branch]
        leaves = [s for s in d.all_structures if s.is_leaf]

        assert len(branches) == 2
        assert len(leaves) == 3

    def test_position_criterion(self):

        def position(structure, index=None, value=None):
            return (np.any(structure.indices()[0] == 6) or
                    np.any(structure.indices()[0] == 8))

        d = Dendrogram.compute(self.data, is_independent=position)

        branches = [s for s in d.all_structures if s.is_branch]
        leaves = [s for s in d.all_structures if s.is_leaf]

        assert len(branches) == 1
        assert len(leaves) == 2

        # Check that leaf that used to contain pixels 1, 2, and 3 is now just
        # part of the main branch.
        assert np.all(branches[0].indices(subtree=False) == np.array([0, 1, 2, 3, 4, 7, 9]))


# Try and reproduce the benchmark tests using a function instead of arguments


from .build_benchmark import BENCHMARKS


@pytest.mark.parametrize(('filename'), BENCHMARKS.keys())
def test_benchmark(filename):

    from astropy.io import fits
    import os

    path = os.path.join(os.path.dirname(__file__),
                        'benchmark_data', filename)

    p = BENCHMARKS[filename]
    data = fits.getdata(path, 1)

    # Now define a function that will test the criteria
    def is_independent_test(structure, index=None, value=None):
        if value is None:
            value = structure.vmin
        if 'min_delta' in p:
            if structure.vmax - value < p['min_delta']:
                return False
        if 'min_npix' in p:
            if len(structure.values()) < p['min_npix']:
                return False
        return True

    d1 = Dendrogram.compute(data,
                            is_independent=is_independent_test,
                            min_value=p['min_value'] if 'min_value' in p else "min")
    d2 = Dendrogram.load_from(path)

    assert d1 == d2

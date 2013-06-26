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

" Test import and export of dendrograms "

import os

import numpy as np

from .. import Dendrogram
from ..structure import Structure
from .test_index import assert_permuted_fancyindex

class TestIO(object):

    def setup_method(self, method):
        n = np.nan
        self.data = np.array([[[n, n, n, n, n, n, n, n],
                               [n, 4, n, n, n, n, n, n],
                               [n, n, n, 1, n, n, 0, 5],
                               [3, n, n, 2, 3, 2, 0, n]],
                              [[n, n, n, n, n, n, n, n],
                               [1, n, n, n, n, n, n, n],
                               [1, n, 1, 1, 0, n, 0, 1],
                               [2, n, n, 1, 3, 1, n, 1]],
                              [[n, 2, 3, 4, n, n, n, n],
                               [1, 1, n, n, n, n, n, n],
                               [n, n, n, n, n, n, 0, 1],
                               [n, n, n, 1, 0, 1, 0, n]]])
        self.test_filename = None

    def teardown_method(self, method):
        if self.test_filename and os.path.exists(self.test_filename):
            os.remove(self.test_filename)

    def compare_dendrograms(self, d1, d2):
        " Helper method that ensures d1 and d2 are equivalent "
        # Do we get the same number of structures?
        assert len(d1.structures_dict) == len(d2.structures_dict)
        # Do we recover the data exactly?
        np.testing.assert_array_equal(d1.data, d2.data)
        # Now check that the structures are the same:
        for idx in d2.structures_dict:
            structure1, structure2 = d1.structures_dict[idx], d2.structures_dict[idx]
            assert_permuted_fancyindex(structure1.indices(), structure2.indices())
            assert np.all(np.sort(structure1.values()) == np.sort(structure2.values()))
            assert type(structure1) == type(structure2)
            # Compare the coordinates and data values of all peak pixels:
            assert structure1.get_peak(subtree=True) == structure2.get_peak(subtree=True)

    # Below are the actual tests for each import/export format:

    def test_hdf5(self):
        self.test_filename = 'astrodendro-test.hdf5'
        d1 = Dendrogram.compute(self.data, verbose=False)
        d1.save_to(self.test_filename)
        d2 = Dendrogram.load_from(self.test_filename)
        self.compare_dendrograms(d1, d2)

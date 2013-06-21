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

import timeit
import os

from ._testdata import data

from .. import Dendrogram


def benchmark_compute():
    print("Data loaded. Starting dendrogram computations...")

    def test1():
        Dendrogram.compute(data, min_npix=4, min_value=1.4, min_delta=0.3)
        print("  Completed an iteration of test1.")

    def test2():
        Dendrogram.compute(data, min_npix=8, min_value=1.4)
        print("  Completed an iteration of test2.")

    num = 3

    t1 = timeit.timeit(test1, number=num) / num
    t2 = timeit.timeit(test2, number=num) / num

    print("test1 average over {num} computations: {result:.3} s".format(num=num, result=t1))
    print("test2 average over {num} computations: {result:.3} s".format(num=num, result=t2))

    print("Total average compute time: {0:.3}s".format((t1 + t2) / 2))


def benchmark_hdf5():
    print("\nGenerating complex dendrogram for HDF5 import/export...")
    d = Dendrogram.compute(data, min_npix=2, min_value=1.4, min_delta=0.01)
    print("Dendrogram generated. Testing import and export...")

    filename = '.astrodendro-hdf5-benchmark.hdf5'
    if os.path.exists(filename):
        os.remove(filename)
    num = 2

    def testHDF5():
        print('Exporting...')
        d.save_to(filename)
        print('Importing...')
        d2 = Dendrogram.load_from(filename)
        os.remove(filename)

    t = timeit.timeit(testHDF5, number=num) / num

    print("Total average export+import time: {0:.3}s".format(t))


if __name__ == '__main__':
    try:
        benchmark_compute()
        benchmark_hdf5()
    except KeyboardInterrupt:
        print("Cancelled.")

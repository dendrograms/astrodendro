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

import numpy as np

from .util import parse_dendrogram


def dendro_export_hdf5(d, filename):
    """Export the dendrogram 'd' to the HDF5 file 'filename'"""
    import h5py
    f = h5py.File(filename, 'w')

    f.attrs['n_dim'] = d.n_dim

    f.create_dataset('newick', data=d.to_newick())

    ds = f.create_dataset('index_map', data=d.index_map, compression=True)
    ds.attrs['CLASS'] = 'IMAGE'
    ds.attrs['IMAGE_VERSION'] = '1.2'
    ds.attrs['IMAGE_MINMAXRANGE'] = [d.index_map.min(), d.index_map.max()]

    ds = f.create_dataset('data', data=d.data, compression=True)
    ds.attrs['CLASS'] = 'IMAGE'
    ds.attrs['IMAGE_VERSION'] = '1.2'
    ds.attrs['IMAGE_MINMAXRANGE'] = [d.data.min(), d.data.max()]

    f.close()


def dendro_import_hdf5(filename):
    """Import 'filename' and construct a dendrogram from it"""
    import h5py

    with h5py.File(filename, 'r') as h5f:
        newick = h5f['newick'].value
        data = h5f['data'].value
        index_map = h5f['index_map'].value

    return parse_dendrogram(newick, data, index_map)

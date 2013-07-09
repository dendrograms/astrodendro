# Licensed under an MIT open source license - see LICENSE

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

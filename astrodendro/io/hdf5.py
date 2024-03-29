# Licensed under an MIT open source license - see LICENSE

import os

from astropy import log

from .util import parse_dendrogram
from .handler import IOHandler

HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'


def is_hdf5(filename, mode='r'):
    if mode == 'r' and os.path.exists(filename):
        fileobj = open(filename, 'rb')
        sig = fileobj.read(8)
        return sig == HDF5_SIGNATURE
    elif filename.lower().endswith(('.hdf5', '.h5')):
        return True
    else:
        return False


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

    for key in d.params.keys():
        f.attrs[key] = d.params[key]

    try:
        f.create_dataset('wcs_header', data=d.wcs.to_header_string())
    except AttributeError:
        pass

    f.close()


def dendro_import_hdf5(filename):
    """Import 'filename' and construct a dendrogram from it"""
    import h5py
    from astropy.wcs.wcs import WCS

    log.debug('Loading HDF5 file from disk...')
    with h5py.File(filename, 'r') as h5f:
        newick = h5f['newick'][()].decode("utf-8")  # str
        data = h5f['data'][:]  # numpy array
        index_map = h5f['index_map'][:]  # numpy array
        params = {}
        if 'min_value' in h5f.attrs:
            params['min_value'] = h5f.attrs['min_value']
            params['min_delta'] = h5f.attrs['min_delta']
            params['min_npix'] = h5f.attrs['min_npix']

        try:
            wcs = WCS(h5f['wcs_header'][()])
        except KeyError:
            wcs = None

    log.debug('Parsing dendrogram...')
    return parse_dendrogram(newick, data, index_map, params, wcs)


HDF5Handler = IOHandler(identify=is_hdf5,
                        export_dendro=dendro_export_hdf5,
                        import_dendro=dendro_import_hdf5)

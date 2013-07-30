# Licensed under an MIT open source license - see LICENSE

from .fits import FITSHandler
from .hdf5 import HDF5Handler


IO_FORMATS = {
    'fits': FITSHandler,
    'hdf5': HDF5Handler
}


def load_dendrogram(filename, format=None):
    if format is not None:
        return IO_FORMATS[format].import_dendro(filename)
    else:
        for io_format in IO_FORMATS:
            if IO_FORMATS[io_format].identify(filename, mode='r'):
                return IO_FORMATS[io_format].import_dendro(filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are 'fits' and 'hdf5')")


def save_dendrogram(dendrogram, filename, format=None):
    if format is not None:
        return IO_FORMATS[format].export_dendro(dendrogram, filename)
    else:
        for io_format in IO_FORMATS:
            if IO_FORMATS[io_format].identify(filename, mode='w'):
                return IO_FORMATS[io_format].export_dendro(dendrogram, filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are 'fits' and 'hdf5')")

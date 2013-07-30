# Licensed under an MIT open source license - see LICENSE

from .fits import FITSHandler
from .hdf5 import HDF5Handler


IO_FORMATS = {
    'fits': FITSHandler,
    'hdf5': HDF5Handler
}


def _valid(formats):
    return " and ".join(["'{0}'".format(key) for key in formats])


def load_dendrogram(filename, format=None):
    if format is not None:
        return IO_FORMATS[format].import_dendro(filename)
    else:
        for handler in IO_FORMATS.values():
            if handler.identify(filename, mode='r'):
                return handler.import_dendro(filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are {0})".format(_valid(IO_FORMATS)))


def save_dendrogram(dendrogram, filename, format=None):
    if format is not None:
        return IO_FORMATS[format].export_dendro(dendrogram, filename)
    else:
        for handler in IO_FORMATS.values():
            if handler.identify(filename, mode='w'):
                return handler.export_dendro(dendrogram, filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are {0})".format(_valid(IO_FORMATS)))

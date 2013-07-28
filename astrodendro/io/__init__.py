# Licensed under an MIT open source license - see LICENSE

from .fits import is_fits, dendro_export_fits, dendro_import_fits
from .hdf5 import is_hdf5, dendro_export_hdf5, dendro_import_hdf5


IO_FORMATS = {
    # name: (identify, export_function, import_function)
    'fits': {'identify': is_fits,
             'export_dendro': dendro_export_fits,
             'import_dendro': dendro_import_fits,
             },
    'hdf5': {'identify': is_hdf5,
             'export_dendro': dendro_export_hdf5,
             'import_dendro': dendro_import_hdf5,
             },
}


def load_dendrogram(filename, format=None):
    for io_format in IO_FORMATS:
        if IO_FORMATS[io_format]['identify'](filename, mode='r'):
            return IO_FORMATS[io_format]['import_dendro'](filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are 'fits' and 'hdf5')")


def save_dendrogram(dendrogram, filename, format=None):
    for io_format in IO_FORMATS:
        if IO_FORMATS[io_format]['identify'](filename, mode='w'):
            return IO_FORMATS[io_format]['export_dendro'](dendrogram, filename)
    raise IOError("Could not automatically identify file format - use the "
                  "format= option to specify which format to use (valid "
                  "options are 'fits' and 'hdf5')")

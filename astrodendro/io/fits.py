# Licensed under an MIT open source license - see LICENSE

import numpy as np

# Import and export


def dendro_export_fits(d, filename):
    """Export the dendrogram 'd' to the FITS file 'filename'"""
    import pyfits
    raise NotImplementedError("FITS export has not yet been implemented.")


def dendro_import_fits(filename):
    """Import 'filename' and construct a dendrogram from it"""
    import pyfits
    from ..dendrogram import Dendrogram
    from ..structure import Structure
    raise NotImplementedError("FITS import has not yet been implemented.")

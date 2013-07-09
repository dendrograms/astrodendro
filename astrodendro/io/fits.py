# Licensed under an MIT open source license - see LICENSE

import numpy as np

from .util import parse_dendrogram
# Import and export


def dendro_export_fits(d, filename):
    """Export the dendrogram 'd' to the FITS file 'filename'"""
    from astropy.io import fits

    hdus = [fits.PrimaryHDU(),
            fits.ImageHDU(d.data),
            fits.ImageHDU(d.index_map),
            fits.ImageHDU(np.array([ord(x) for x in  d.to_newick()]))]
    hdulist = fits.HDUList(hdus)

    hdulist.writeto(filename, clobber=True)


def dendro_import_fits(filename):
    """Import 'filename' and construct a dendrogram from it"""
    from astropy.io import fits

    with fits.open(filename) as hdus:
        data = hdus[1].data
        index_map = hdus[2].data
        newick = ''.join(chr(x) for x in hdus[3].data.flat)

    return parse_dendrogram(newick, data, index_map)

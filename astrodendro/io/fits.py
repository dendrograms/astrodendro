# Licensed under an MIT open source license - see LICENSE

import os

import numpy as np

from .util import parse_dendrogram
from .handler import IOHandler

# Import and export

# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")


def is_fits(filename, mode='r'):
    if mode == 'r' and os.path.exists(filename):
        fileobj = open(filename, 'rb')
        sig = fileobj.read(30)
        return sig == FITS_SIGNATURE
    elif filename.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz')):
        return True
    else:
        return False


def dendro_export_fits(d, filename):
    """Export the dendrogram 'd' to the FITS file 'filename'"""
    from astropy.io import fits

    try:
        primary_hdu = fits.PrimaryHDU(header=d.wcs.to_header())
    except AttributeError:
        primary_hdu = fits.PrimaryHDU()

    hdus = [primary_hdu,
            fits.ImageHDU(d.data),
            fits.ImageHDU(d.index_map),
            fits.ImageHDU(np.array([ord(x) for x in d.to_newick()]))]

    hdulist = fits.HDUList(hdus)

    hdulist.writeto(filename, clobber=True)


def dendro_import_fits(filename):
    """Import 'filename' and construct a dendrogram from it"""
    from astropy.io import fits
    from astropy.wcs.wcs import WCS

    with fits.open(filename) as hdus:
        try:
            wcs = WCS(hdus[0].header)
        except AttributeError:
            wcs = None
        data = hdus[1].data
        index_map = hdus[2].data
        newick = ''.join(chr(x) for x in hdus[3].data.flat)

    return parse_dendrogram(newick, data, index_map, wcs)


FITSHandler = IOHandler(identify=is_fits,
                        export_dendro=dendro_export_fits,
                        import_dendro=dendro_import_fits)

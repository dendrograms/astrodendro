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

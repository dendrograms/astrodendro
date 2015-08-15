from __future__ import print_function, division

import os
from ..fits import dendro_import_fits

DATA = os.path.join(os.path.dirname(__file__), 'data')


def test_import_old():
    # Check that we are backward-compatible
    dendro_import_fits(os.path.join(DATA, 'dendro_old.fits'))


def test_import():
    dendro_import_fits(os.path.join(DATA, 'dendro.fits'))

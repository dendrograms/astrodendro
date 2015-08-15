from __future__ import print_function, division

import os
from ..hdf5 import dendro_import_hdf5

DATA = os.path.join(os.path.dirname(__file__), 'data')


def test_import():
    dendro_import_hdf5(os.path.join(DATA, 'dendro.hdf5'))

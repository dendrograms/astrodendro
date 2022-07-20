# Licensed under an MIT open source license - see LICENSE

" Test import and export of dendrograms "

import pytest
import numpy as np

from .. import Dendrogram
from .test_index import assert_permuted_fancyindex

from astropy.wcs import WCS


class TestIO(object):

    def setup_method(self, method):
        n = np.nan
        self.data = np.array([[[n, n, n, n, n, n, n, n],
                               [n, 4, n, n, n, n, n, n],
                               [n, n, n, 1, n, n, 0, 5],
                               [3, n, n, 2, 3, 2, 0, n]],
                              [[n, n, n, n, n, n, n, n],
                               [1, n, n, n, n, n, n, n],
                               [1, n, 1, 1, 0, n, 0, 1],
                               [2, n, n, 1, 3, 1, n, 1]],
                              [[n, 2, 3, 4, n, n, n, n],
                               [1, 1, n, n, n, n, n, n],
                               [n, n, n, n, n, n, 0, 1],
                               [n, n, n, 1, 0, 1, 0, n]]])

    def compare_dendrograms(self, d1, d2):
        " Helper method that ensures d1 and d2 are equivalent "
        # Do we get the same number of structures?
        assert len(d1) == len(d2)
        # Do we recover the data exactly?
        np.testing.assert_array_equal(d1.data, d2.data)
        # Now check that the structures are the same:
        for s in d2:
            idx = s.idx
            structure1, structure2 = d1[idx], d2[idx]
            assert_permuted_fancyindex(structure1.indices(subtree=False),
                                       structure2.indices(subtree=False))
            assert np.all(np.sort(structure1.values(subtree=False)) ==
                          np.sort(structure2.values(subtree=False)))
            assert isinstance(structure1, type(structure2))
            # Compare the coordinates and data values of all peak pixels:
            assert structure1.get_peak(subtree=True) == \
                structure2.get_peak(subtree=True)

            assert structure2._tree_index is not None

    # Below are the actual tests for each import/export format:

    def test_hdf5(self, tmp_path):
        test_filename = tmp_path / 'astrodendro-test.hdf5'
        d1 = Dendrogram.compute(self.data, verbose=False)
        d1.save_to(test_filename, format='hdf5')
        d2 = Dendrogram.load_from(test_filename, format='hdf5')
        self.compare_dendrograms(d1, d2)

    def test_fits(self, tmp_path):
        test_filename = tmp_path / 'astrodendro-test.fits'
        d1 = Dendrogram.compute(self.data, verbose=False)
        d1.save_to(test_filename, format='fits')
        d2 = Dendrogram.load_from(test_filename, format='fits')
        self.compare_dendrograms(d1, d2)

    def test_hdf5_auto(self, tmp_path):

        test_filename = tmp_path / 'astrodendro-test.hdf5'
        test_filename_noext = tmp_path / 'astrodendro-test'

        d1 = Dendrogram.compute(self.data, verbose=False)

        # recognize from extension
        d1.save_to(test_filename)

        # no way to tell
        with pytest.raises(IOError):
            d1.save_to(test_filename_noext)

        # no way to tell, so have to explicitly give format
        d1.save_to(test_filename_noext, format='hdf5')

        # recognize from extension
        Dendrogram.load_from(test_filename)

        # recognize from signature
        Dendrogram.load_from(test_filename_noext)

    def test_fits_auto(self, tmp_path):

        test_filename = tmp_path / 'astrodendro-test.fits'
        test_filename_noext = tmp_path / 'astrodendro-test'

        d1 = Dendrogram.compute(self.data, verbose=False)

        # recognize from extension
        d1.save_to(test_filename)

        # no way to tell
        with pytest.raises(IOError):
            d1.save_to(test_filename_noext)

        # no way to tell, so have to explicitly give format
        d1.save_to(test_filename_noext, format='fits')

        # recognize from extension
        Dendrogram.load_from(test_filename)

        # recognize from signature
        Dendrogram.load_from(test_filename_noext)

    def test_hdf5_with_wcs(self, tmp_path):
        test_filename = tmp_path / 'astrodendro-test-wcs.hdf5'
        test_wcs = WCS(header=dict(cdelt1=1, crval1=0, crpix1=1,
                                   cdelt2=2, crval2=0, crpix2=1,
                                   cdelt3=3, crval3=0, crpix3=1))

        d1 = Dendrogram.compute(self.data, verbose=False, wcs=test_wcs)
        d1.save_to(test_filename, format='hdf5')
        d2 = Dendrogram.load_from(test_filename, format='hdf5')

        assert d2.wcs.to_header_string() == d1.wcs.to_header_string()

    def test_fits_with_wcs(self, tmp_path):
        test_filename = tmp_path / 'astrodendro-test-wcs.fits'
        test_wcs = WCS(header=dict(cdelt1=1, crval1=0, crpix1=1,
                                   cdelt2=2, crval2=0, crpix2=1,
                                   cdelt3=3, crval3=0, crpix3=1))
        d1 = Dendrogram.compute(self.data, verbose=False, wcs=test_wcs)
        d1.save_to(test_filename, format='fits')
        d2 = Dendrogram.load_from(test_filename, format='fits')

        assert d2.wcs.to_header_string() == d1.wcs.to_header_string()

    @pytest.mark.parametrize('ext', ('fits', 'hdf5'))
    def test_reload_retains_dendro_reference(self, ext, tmp_path):
        # regression test for issue 106

        d1 = Dendrogram.compute(self.data, verbose=False)

        test_filename = tmp_path / f'astrodendro-test.{ext}'

        d1.save_to(test_filename)
        d2 = Dendrogram.load_from(test_filename)

        for s in d1:
            np.testing.assert_array_equal(d2[s.idx].get_mask(subtree=True),
                                          d1[s.idx].get_mask(subtree=True))

import pytest
import numpy as np

from astropy import units as u

from ..flux import compute_flux


COMBINATIONS = \
    [
        (np.array([1, 2, 3]) * u.Jy, u.Jy, {}, 6 * u.Jy),
        (np.array([1, 2, 3]) * u.mJy, u.Jy, {}, 0.006 * u.Jy),
        (np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.Hz, u.Jy, {}, 6e23 * u.Jy),
        (np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.micron, u.Jy, {'wavelength': 2 * u.micron}, 8005538284.75565 * u.Jy),
        (np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy, {'spatial_scale': 2 * u.arcsec}, 24. * u.Jy),
        (np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, {'spatial_scale': 2 * u.arcsec, 'beam_major': 1 * u.arcsec, 'beam_minor': 0.5 * u.arcsec}, 42.36166269526079 * u.Jy),
        (np.array([1, 2, 3]) * u.K, u.Jy, {'spatial_scale': 2 * u.arcsec, 'beam_major': 1 * u.arcsec, 'beam_minor': 0.5 * u.arcsec, 'wavelength': 2 * u.mm}, 0.38941636582186634 * u.Jy),
        (np.array([1, 2, 3]) * u.K, u.Jy, {'spatial_scale': 2 * u.arcsec, 'beam_major': 1 * u.arcsec, 'beam_minor': 0.5 * u.arcsec, 'wavelength': 100 * u.GHz}, 0.17331365650395836 * u.Jy),
    ]


@pytest.mark.parametrize(('input_quantities', 'output_unit', 'keywords', 'output'), COMBINATIONS)
def test_compute_flux(input_quantities, output_unit, keywords, output):
    q = compute_flux(input_quantities, output_unit, **keywords)
    np.testing.assert_allclose(q.value, output.value, rtol=1e-6)
    assert q.unit == output.unit


def test_monochromatic_wav_missing():
    with pytest.raises(ValueError, match='wavelength is needed to convert from erg'):
        compute_flux(np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.micron, u.Jy)


def test_monochromatic_wav_invalid_units():
    with pytest.raises(ValueError, match='wavelength should be a physical length'):
        compute_flux(np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.micron, u.Jy, wavelength=3 * u.L)


def test_surface_brightness_scale_missing():
    with pytest.raises(ValueError, match='spatial_scale is needed to convert from Jy'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy)


def test_surface_brightness_invalid_units():
    with pytest.raises(ValueError, match='spatial_scale should be an angle'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy, spatial_scale=3 * u.m)


def test_per_beam_scale_missing():

    with pytest.raises(ValueError, match='spatial_scale is needed to convert from Jy / beam to Jy'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, beam_major=3 * u.arcsec, beam_minor=2. * u.arcsec)

    with pytest.raises(ValueError, match='beam_major is needed to convert from Jy / beam to Jy'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, spatial_scale=3 * u.arcsec, beam_minor=2. * u.arcsec)

    with pytest.raises(ValueError, match='beam_minor is needed to convert from Jy / beam to Jy'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, spatial_scale=3 * u.arcsec, beam_major=2. * u.arcsec)


def test_per_beam_invalid_units():

    with pytest.raises(ValueError, match='beam_major should be an angle'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, spatial_scale=3 * u.arcsec, beam_major=3 * u.m, beam_minor=2. * u.arcsec)

    with pytest.raises(ValueError, match='beam_minor should be an angle'):
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.beam, u.Jy, spatial_scale=3 * u.arcsec, beam_major=3 * u.arcsec, beam_minor=2. * u.m)

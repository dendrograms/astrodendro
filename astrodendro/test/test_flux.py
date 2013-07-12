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
(np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy, {'pixel_scale': 2 * u.arcsec}, 24. * u.Jy),
    ]


@pytest.mark.parametrize(('input_quantities', 'output_unit', 'keywords', 'output'), COMBINATIONS)
def test_compute_flux(input_quantities, output_unit, keywords, output):
    q = compute_flux(input_quantities, output_unit, **keywords)
    np.testing.assert_allclose(q.value, output.value)
    assert q.unit == output.unit


def test_monochromatic_wav_missing():
    with pytest.raises(ValueError) as exc:
        compute_flux(np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.micron, u.Jy)
    assert exc.value.args[0] == 'Wavelength is needed to convert from erg / (cm2 micron s) to Jy'


def test_monochromatic_wav_invalid_units():
    with pytest.raises(ValueError) as exc:
        compute_flux(np.array([1, 2, 3]) * u.erg / u.cm ** 2 / u.s / u.micron, u.Jy, wavelength=3 * u.L)
    assert exc.value.args[0] == 'Wavelength should be a physical length'


def test_surface_brightness_scale_missing():
    with pytest.raises(ValueError) as exc:
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy)
    assert exc.value.args[0] == 'Pixel scale is needed to convert from Jy / arcsec2 to Jy'


def test_surface_brightness_scale_missing():
    with pytest.raises(ValueError) as exc:
        compute_flux(np.array([1, 2, 3]) * u.Jy / u.arcsec ** 2, u.Jy, pixel_scale=3 * u.m)
    assert exc.value.args[0] == 'Pixel scale should be an angle'

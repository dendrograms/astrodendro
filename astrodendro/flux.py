import warnings

import numpy as np

from astropy import units as u
from astropy.constants import si

class UnitMetadataWarning(UserWarning):
   pass

def quantity_sum(quantities):
    """
    In Astropy 0.3, np.sum will do the right thing for quantities, but in the mean time we need a workaround.
    """
    return np.sum(quantities.value) * quantities.unit


def compute_flux(input_quantities, output_unit, wavelength=None, spatial_scale=None,
                 velocity_scale=None, beam_major=None, beam_minor=None):
    """
    Given a set of flux values in arbitrary units, find the total flux in a
    specific set of units.

    Parameters
    ----------
    input_quantities : :class:`~astropy.units.quantity.Quantity` instance
        A `astropy.units.quantity.Quantity` instance containing an array of
        flux values to be summed.
    output_unit : :class:`~astropy.units.core.Unit` instance
        The final unit to give the total flux in (should be equivalent to Jy)
    wavelength : :class:`~astropy.units.quantity.Quantity` instance
        The wavelength of the data (required if converting e.g.
        ergs/cm^2/s/micron to Jy)
    spatial_scale : :class:`~astropy.units.quantity.Quantity` instance
        The pixel scale of the data (should be an angle)
    velocity_scale : :class:`~astropy.units.quantity.Quantity` instance
        The pixel scale of the data (should be a velocity)
    beam_major : :class:`~astropy.units.quantity.Quantity` instance
        The beam major full width at half_maximum (FWHM)
    beam_minor : :class:`~astropy.units.quantity.Quantity` instance
        The beam minor full width at half_maximum (FWHM)
    """

    # Start off by finding the total flux in Jy

    if input_quantities.unit.is_equivalent(u.Jy):  # Fnu

        # Simply sum up the values and convert to output unit
        total_flux = quantity_sum(input_quantities).to(u.Jy)

    elif input_quantities.unit.is_equivalent(u.erg / u.cm ** 2 / u.s / u.m):  # Flambda

        if wavelength is not None and not wavelength.unit.is_equivalent(u.m):
            raise ValueError("wavelength should be a physical length")

        # Find the frequency
        if wavelength is None:
            raise ValueError("wavelength is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find frequency
        nu = si.c / wavelength

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * wavelength / nu).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.MJy / u.sr):  # surface brightness (Fnu)

        if spatial_scale is not None and not spatial_scale.unit.is_equivalent(u.degree):
            raise ValueError("spatial_scale should be an angle")

        if spatial_scale is None:
            raise ValueError("spatial_scale is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find the area of a pixel as a solid angle
        pixel_area = (spatial_scale ** 2)

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * pixel_area).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.Jy / u.beam):

        if spatial_scale is not None and not spatial_scale.unit.is_equivalent(u.degree):
            raise ValueError("spatial_scale should be an angle")

        if spatial_scale is None:
            raise ValueError("spatial_scale is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_major is not None and not beam_major.unit.is_equivalent(u.degree):
            raise ValueError("beam_major should be an angle")

        if beam_major is None:
            raise ValueError("beam_major is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_minor is not None and not beam_minor.unit.is_equivalent(u.degree):
            raise ValueError("beam_minor should be an angle")

        if beam_minor is None:
            raise ValueError("beam_minor is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find the beam area
        beams_per_pixel = spatial_scale ** 2 / (beam_minor * beam_major * 1.1331) * u.beam

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * beams_per_pixel).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.K):

        if spatial_scale is not None and not spatial_scale.unit.is_equivalent(u.degree):
            raise ValueError("spatial_scale should be an angle")

        if spatial_scale is None:
            raise ValueError("spatial_scale is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_major is not None and not beam_major.unit.is_equivalent(u.degree):
            raise ValueError("beam_major should be an angle")

        if beam_major is None:
            raise ValueError("beam_major is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_minor is not None and not beam_minor.unit.is_equivalent(u.degree):
            raise ValueError("beam_minor should be an angle")

        if beam_minor is None:
            raise ValueError("beam_minor is needed to convert from {0} to Jy".format(input_quantities.unit))

        if wavelength is not None and not wavelength.unit.is_equivalent(u.m, equivalencies=u.spectral()):
            raise ValueError("wavelength should be a physical length")

        # Find the frequency
        if wavelength is None:
            raise ValueError("wavelength is needed to convert from {0} to Jy".format(input_quantities.unit))

        warnings.warn("'Kelvin' units interpreted as main beam brightness temperature.",
                      UnitMetadataWarning)

        # Find frequency
        nu = wavelength.to(u.Hz, equivalencies=u.spectral())

        # Angular area of beam. Conversion between 2D Gaussian FWHM and effective area comes from https://github.com/radio-astro-tools/radio_beam/blob/bc906c38a65e85c6a894ee81519a642665e50f7c/radio_beam/beam.py#L8
        omega_beam = np.pi * 2 / (8*np.log(2)) * beam_major * beam_minor

        # Find the beam area
        beams_per_pixel = spatial_scale ** 2 / omega_beam * u.beam

        # Convert input quantity to Fnu in Jy
        # Implicitly, this equivalency gives the Janskys in a single beam, so we make this explicit by dividing out a beam
        jansky_per_beam = input_quantities.to(u.Jy, 
            equivalencies=u.brightness_temperature(omega_beam, nu)) / u.beam

        q = jansky_per_beam * beams_per_pixel

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    else:

        raise ValueError("Flux units {0} not yet supported".format(input_quantities.unit))

    if not output_unit.is_equivalent(u.Jy):
        raise ValueError("output_unit has to be equivalent to Jy")
    else:
        return total_flux.to(output_unit)

import numpy as np

from astropy import units as u
from astropy.constants import si


def quantity_sum(quantities):
    """
    In Astropy 0.3, np.sum will do the right thing for quantities, but in the mean time we need a workaround
    """
    return np.sum(quantities.value) * quantities.unit


def compute_flux(input_quantities, output_unit, wavelength=None, pixel_scale=None,
                 velocity_scale=None, beam_major=None, beam_minor=None):
    """
    Given a set of flux values in arbitrary units, find the total flux in a
    specific set of units.

    Parameters
    ----------
    input_quantities : `~astropy.units.quantity.Quantity` instance
        A `~astropy.units.quantity.Quantity` instance containing an array of
        flux values to be summed.
    output_unit : `~astropy.units.core.Unit` instance
        The final unit to give the total flux in (should be equivalent to Jy)
    wavelength : `~astropy.units.quantity.Quantity` instance
        The wavelength of the data (required if converting e.g.
        ergs/cm^2/s/micron to Jy)
    pixel_scale : `~astropy.units.quantity.Quantity` instance
        The pixel scale of the data (should be an angle)
    velocity_scale : `~astropy.units.quantity.Quantity` instance
        The pixel scale of the data (should be a velocity)
    beam_major : `~astropy.units.quantity.Quantity` instance
        The beam major full width at half_maximum (FWHM)
    beam_minor : `~astropy.units.quantity.Quantity` instance
        The beam minor full width at half_maximum (FWHM)
    """

    if pixel_scale is not None and not pixel_scale.unit.is_equivalent(u.degree):
        raise ValueError("Pixel scale should be an angle")

    if velocity_scale is not None and not velocity_scale.unit.is_equivalent(u.m / u.s):
        raise ValueError("Pixel scale should be an angle")

    if wavelength is not None and not wavelength.unit.is_equivalent(u.m):
        raise ValueError("Wavelength should be a physical length")

    if beam_major is not None and not beam_major.unit.is_equivalent(u.degree):
        raise ValueError("Beam major FWHM should be an angle")

    if beam_minor is not None and not beam_minor.unit.is_equivalent(u.degree):
        raise ValueError("Beam minor FWHM should be an angle")

    # Start off by finding the total flux in Jy

    if input_quantities.unit.is_equivalent(u.Jy):  # Fnu

        # Simply sum up the values and convert to output unit
        total_flux = quantity_sum(input_quantities).to(u.Jy)

    elif input_quantities.unit.is_equivalent(u.erg / u.cm ** 2 / u.s / u.m):  # Flambda

        # Find the frequency
        if wavelength is None:
            raise ValueError("Wavelength is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find frequency
        nu = si.c / wavelength

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * wavelength / nu).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.MJy / u.sr):  # surface brightness (Fnu)

        if pixel_scale is None:
            raise ValueError("Pixel scale is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find the area of a pixel as a solid angle
        pixel_area = (pixel_scale ** 2)

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * pixel_area).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.Jy / u.beam):

        if pixel_scale is None:
            raise ValueError("Pixel scale is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_major is None:
            raise ValueError("Beam major FWHM is needed to convert from {0} to Jy".format(input_quantities.unit))

        if beam_minor is None:
            raise ValueError("Beam minor FWHM is needed to convert from {0} to Jy".format(input_quantities.unit))

        # Find the beam area
        beams_per_pixel = (beam_minor * beam_major * 1.1331 / pixel_scale ** 2) * u.beam

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * beams_per_pixel).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    else:

        raise ValueError("Flux units {0} not yet supported".format(input_quantities.unit))

    if not output_unit.is_equivalent(u.Jy):
        raise ValueError("Output unit has to be equivalent to Jy")
    else:
        return total_flux.to(output_unit)

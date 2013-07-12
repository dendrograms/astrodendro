import numpy as np

from astropy import units as u
from astropy.constants import si


def quantity_sum(quantities):
    """
    In Astropy 0.3, np.sum will do the right thing for quantities, but in the mean time we need a workaround
    """
    return np.sum(quantities.value) * quantities.unit


def compute_flux(input_quantities, output_unit, wavelength=None, pixel_scale=None, velocity_scale=None):
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
    """

    # Start off by finding the total flux in Jy

    if input_quantities.unit.is_equivalent(u.Jy):  # monochromatic (frequency)

        # Simply sum up the values and convert to output unit
        total_flux = quantity_sum(input_quantities).to(u.Jy)

    elif input_quantities.unit.is_equivalent(u.erg / u.cm ** 2 / u.s / u.m):

        # Find the frequency
        if wavelength is None:
            raise ValueError("Wavelength is needed to convert from {0} to Jy".format(input_quantities.unit))
        elif not wavelength.unit.is_equivalent(u.m):
            raise ValueError("Wavelength should be a physical length")

        # Find frequency
        nu = si.c / wavelength

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * wavelength / nu).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    elif input_quantities.unit.is_equivalent(u.MJy / u.sr):

        # Find the frequency
        if pixel_scale is None:
            raise ValueError("Pixel scale is needed to convert from {0} to Jy".format(input_quantities.unit))
        elif not pixel_scale.unit.is_equivalent(u.degree):
            raise ValueError("Pixel scale should be an angle")

        # Find the area of a pixel as a solid angle
        pixel_area = (pixel_scale ** 2)

        # Convert input quantity to Fnu in Jy
        q = (input_quantities * pixel_area).to(u.Jy)

        # Find total flux in Jy
        total_flux = quantity_sum(q)

    else:

        raise ValueError("Units {0} not yet supported".format(input_quantities.unit))

    if not output_unit.is_equivalent(u.Jy):
        raise ValueError("Output unit has to be equivalent to Jy")
    else:
        return total_flux.to(output_unit)

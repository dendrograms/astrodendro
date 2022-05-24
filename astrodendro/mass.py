import numpy as np

from astropy import units as u
from astropy.constants import si


def quantity_sum(quantities):
    """
    In Astropy 0.3, np.sum will do the right thing for quantities, but in the mean time we need a workaround.
    """
    return np.sum(quantities.value) * quantities.unit


def compute_mass(input_quantities, output_unit, spatial_scale=None):
    """
    Given a set of density values in Msun/pc**3 or g/cm**3 units, 
    find the total mass in a specific set of units.

    Parameters
    ----------
    input_quantities : `~astropy.units.quantity.Quantity` instance
        A `~astropy.units.quantity.Quantity` instance containing an array of
        density values to be summed.
    output_unit : `~astropy.units.core.Unit` instance
        The final unit to give the total mass in (should be equivalent to Msun)
    spatial_scale : `~astropy.units.quantity.Quantity` instance
        The pixel scale of the data
    """

    # Start off by finding the total mass in Msun

    if input_quantities.unit.is_equivalent(u.Msun):

        # Simply sum up the values and convert to output unit
        total_flux = quantity_sum(input_quantities).to(u.Msun)

    elif input_quantities.unit.is_equivalent(u.g / u.cm ** 3):

        # Convert volume pixel to cm**3
        if spatial_scale is None:
            print("WARNING: spatial_scale not recognized, assumed as cm")
            pixel_volume = u.cm ** 3          

        elif spatial_scale.unit.is_equivalent(u.kpc):
            pixel_volume = ((spatial_scale.to(u.cm)) ** 3)

        elif spatial_scale.unit.is_equivalent(u.pc):        
            pixel_volume = ((spatial_scale.to(u.cm)) ** 3)

        elif spatial_scale.unit.is_equivalent(u.cm):
            pixel_volume = (spatial_scale ** 3)

        else:
            print("WARNING: spatial_scale not recognized, assumed cm")
            pixel_volume = u.cm ** 3            

        # Find total mass in Msun
        total_mass = quantity_sum(input_quantities) * pixel_volume# / (si.M_sun * 1e3)


    elif input_quantities.unit.is_equivalent(u.Msun / u.pc ** 3):

        # Convert volume pixel to pc**3
        if spatial_scale is None:
            print("WARNING: spatial_scale not recognized, assumed pc")
            pixel_volume = u.pc ** 3

        elif spatial_scale.unit.is_equivalent(u.kpc):
            pixel_volume = ((spatial_scale.to(u.pc)) ** 3)

        elif spatial_scale.unit.is_equivalent(u.pc):        
            pixel_volume = (spatial_scale ** 3)

        elif spatial_scale.unit.is_equivalent(u.cm):
            pixel_volume = ((spatial_scale.to(u.pc)) ** 3)

        else:
            print("WARNING: spatial_scale not recognized, assumed pc")
            pixel_volume = u.pc ** 3            

        # Find total mass in Msun
        total_mass = quantity_sum(input_quantities) * pixel_volume

    else:

        print("WARNING: data unit not recognized. Providing direct sum of values.")
        output_unit = input_quantities.unit
        total_mass = quantity_sum(input_quantities)
        return total_mass.to(output_unit)

    if not output_unit.is_equivalent(u.Msun):
        raise ValueError("output_unit has to be equivalent to Msun")
    else:
        return total_mass.to(output_unit)

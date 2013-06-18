"""
Proposed interface for extracting quantities from Structures

We avoid dependencies on the astrodendro API for as long as possible.
Currently, only the last method assumes anything about astrodendro
"""
import numpy as np


class ScalarStatistic(object):
    #This class does all of the heavy computation
    def __init__(self, values, indices):
        """
        Compute pixel-level statistics from
        a scalar field, sampled at specific locations

        Parameters
        ----------
        values : 1D ndarray
                data values to use
        indices: tuple of 1D arrays
                 Location of each element of values.
                 The ith array in the tuple describes the
                 ith positional dimension
        """
        pass

    def mom0(self):
        pass

    def mom1(self):
        pass

    def mom2(self):
        pass

    def mom2_along(self, direction):
        pass

    def count(self):
        pass

    def surface_area(self):
        pass

    def perimeter(self, plane=None):
        pass


class VectorStatistic(object):
    def __init__(self, values_tuple, indices):
        pass

    def divergence(self):
        pass

    def curl(self):
        pass


class PPVStatistic(object):
    #this class doesn't do heavy lifting directly
    #instead, it uses the stat object

    fields = ('mass', 'flux')  # etc

    def __init__(self, stat, metadata):

        """Combine metadata with a ScalarStatistic describing a PPV
        structure extract semantically-meaningful quantities

        Parameters
        ----------
        stat :   ScalarStatistic instance
        metadata : dict
                 Key-value paris of metadata, like 'dx' for
                 the angular size of a pixel
        """
        #init would include a check to make
        #sure metadata has the neeced keys
        #and raises a helpful exception if not
        #
        #furthermore, if metadata values are astropy quantities,
        #this check can make sure everything is in the right units
        pass

    def mass(self):
        dx_linear = np.radians(metadata['dx']) * metadata['dist']
        cell_size = dx_linear ** 2 * dv

        return self.stat.mom0 * cell_size * metadata['xfactor']

    def flux(self):
        pass

    def sky_maj(self):
        pass

    def sky_min(self):
        pass

    def sky_rad(self):
        pass

    def vrms(self):
        pass

    def sky_deconvolved_maj(self):
        pass

    def sky_deconvolved_min(self):
        pass

    def sky_deconvolved_rad(self):
        pass

    def sky_pa(self):
        pass

    def virial(self):
        pass

    def pressure(self):
        pass


class PPPStatistics(object):

    def __init__(self, rhostat, vstat, metadata):
        """
        Combine metadata with density and velocity
        statistics, to compute semantically meaningful
        quantities

        Parameters
        ----------
        rhostat : ScalarStatistic instance
        vstat : VectorStatistic instance
        """

    def mass(self):
        pass

    def volume(self):
        pass

    def surface_area(self):
        pass

    def virial(self):
        pass

    def vrms(self):
        pass

    def vz_rms(self):
        pass

    def pressure_vz(self):
        pass

    def pressure(self):
        pass


def ppv_catalog(structures, metadata, fields=None):
    # this is the first class with a dependency on astrodendro.
    # it assumes:
    #
    # structures have a .values and .indices property
    # if a dendrogram is passed as the first argument,
    #    it has an __iter__ that iterates over Structures.
    #    Other classes could construct Structures, however
    """
    Iterate over a collection of PPV structures,
    extracting several quantities from each, and building
    a catalog

    Parameters
    ----------
    structures : Iterable of Structures
         The structures to catalog. This could
         be a dendrogram, but need not be

    metadata : dict of metadata
    fields : list of strings, optional
             The quantities to extract. If not provided,
             defaults to PPVStatistic.fields
    """
    result = []
    fields = fields or PPVStatistic.fields

    for struct in structures:
        stat = ScalarStatistic(struct.values, struct.indices)
        stat = PPVStatistic(stat, metadata)
        #extract each field by looking up + calling the method
        result.append(dict((lbl, getattr(stat, lbl)())
                           for lbl in fields))

    #return a list of dicts (like this)?
    #or a numpy recarray?
    return result

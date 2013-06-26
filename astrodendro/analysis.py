from types import FunctionType
import warnings

import numpy as np
from astropy.units import Quantity, rad
from astropy.table import Table

from .structure import Structure

__all__ = ['ppv_catalog', 'pp_catalog']


def _qsplit(q):
    """Split a potential astropy Quantity into unit/quantity"""
    if isinstance(1 * q, Quantity):
        return q.unit, q.value

    return 1, q


def _unit(q):
    """Return the units associated with a number, array, unit, or Quantity"""
    if isinstance(1 * q, Quantity):
        return (1 * q).unit


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
        self.values = values.astype(np.float)
        self.indices = indices

    def mom0(self):
        """The sum of the values"""
        return np.nansum(self.values)

    def mom1(self):
        """The intensity-weighted mean position"""
        m0 = self.mom0()
        return [np.nansum(i * self.values) / m0
                for i in self.indices]

    def mom2(self):
        """The intensity-weighted covariance matrix"""
        mom1 = self.mom1()
        mom0 = self.mom0()
        v = self.values / mom0

        nd = len(self.indices)
        zyx = tuple(i - m for i, m in zip(self.indices, mom1))

        result = np.zeros((nd, nd))

        for i in range(nd):
            result[i, i] = np.nansum(v * zyx[i] ** 2)
            for j in range(i + 1, nd):
                result[i, j] = result[j, i] = np.nansum(v * zyx[i] * zyx[j])

        return result

    def mom2_along(self, direction):
        """Intensity-weighted variance/covariance along 1 or more directions

        Parameters
        ----------
        direction : array like
                  One or more set of direction vectors. Need not be normalized

        Returns
        -------
        The variance (or co-variance matrix) of the data along
        the specified direction(s).
        """
        w = np.atleast_2d(direction).astype(np.float)
        for row in w:
            row /= np.linalg.norm(row)

        result = np.dot(np.dot(w, self.mom2()), w.T)
        if result.size == 1:
            result = np.asscalar(result)
        return result

    def paxes(self):
        """
        The principal axes of the data (direction of greatest elongation)

        Returns
        -------
        Ordered list of ndarrays

        Each array is a normalized direction vector. The arrays
        are sorted in decreasing order of elongation of the data
        """
        mom2 = self.mom2()
        w, v = np.linalg.eig(mom2)
        order = np.argsort(w)

        return tuple(v[:, o] for o in order[::-1])

    def projected_paxes(self, axes):
        """
        The principal axes of a projection of the data onto a subspace

        Paramters
        ---------
        axes : array-like, (nnew, nold)
               The projection to take. Each row defines a unit vector in
               the new coordinate system

        Returns
        --------
        tuple of arrays (nnew items)

        The ordered principal axes in the new space
        """
        mom2 = self.mom2_along(axes)
        w, v = np.linalg.eig(mom2)
        order = np.argsort(w)

        return tuple(v[:, o] for o in order[::-1])

    def count(self):
        """Number of elements in the dataset"""
        return self.values.size

    def surface_area(self):
        raise NotImplementedError

    def perimeter(self, plane=None):
        raise NotImplementedError


class VectorStatistic(object):
    def __init__(self, values_tuple, indices):
        raise NotImplementedError

    def divergence(self):
        raise NotImplementedError

    def curl(self):
        raise NotImplementedError


class MetaData(object):
    """A descriptor to wrap around metadata dictionaries

    Let's classes reference self.x instead of self.metadata['x'],
    """
    def __init__(self, key, description, default=1, strict=False):
        """
        Parameters
        ----------
        key : str
               Metadata name.
        description : str
               What the quantity describes
        default : scalar
               Default value if metadata not provided
        strict : bool
               If True, raise KeyError if metadata not provided.
               This overrides default
        """
        self.key = key
        self.description = description or 'no description'
        self.default = default
        self.strict = strict

    def __get__(self, instance, type=None):
        if instance is None:
            return self

        if self.strict and self.key not in instance.metadata:
            raise KeyError("Required metadata item not found: %s" % self)
        return instance.metadata.get(self.key, self.default)

    def __str__(self):
        return "%s (%s)" % (self.key, self.description)


def _missing_metadata(cl, md):
    """Find missing metadata entries in a metadata dict

    Paramters
    ---------
    cls : Class with MetaData descriptors
    md : metadata dictionary
    """
    result = []
    attrs = [getattr(cl, t) for t in dir(cl)]
    return [m for m in attrs if isinstance(m, MetaData)
            and m.key not in md]


def _warn_missing_metadata(cl, md, verbose=True):
    missing = _missing_metadata(cl, md)
    if len(missing) == 0:
        return

    required = [m for m in missing if m.strict]
    if len(required):
        raise RuntimeError(
            "The following missing metadata items are required:\n\t" +
            "\n\t".join(str(m) for m in required))

    if not verbose:
        return

    for m in missing:
        warnings.warn("Missing Metadata:\n\t %s\n\t Defaulting to %s=%s" %
                      (m, m.key, m.default))


class SpatialBase(object):
    dx = MetaData('dx', 'Angular length of a pixel')
    bmaj = MetaData('bmaj', 'Beam major axis, sigma', default=0)
    bmin = MetaData('bmin', 'Beam minor axis, sigma', default=0)
    bunit = MetaData('bunit', 'Unit of intensity')
    dist = MetaData('dist', 'Distance')
    wcs = MetaData('wcs', 'WCS object', default=None)
    wcs_origin = MetaData('wcs_origin', 'origin (1=FITS standard, 0=numpy)',
                          default=1)

    def luminosity(self):
        """Integrated luminosity

        sum(v_i * dx_linear^2 * dv)
        """
        #disambiguate between degree/radian dx
        #if astropy unit is used
        try:
            fac = (1 * self.dx).unit.to(rad)
            fac /= (1 * self.dx).unit
        except AttributeError:
            # metadata not a quantity. Assuming dx=degrees
            fac = np.radians(1)
        return self.dist ** 2 * self.flux() * fac ** 2

    def _sky_paxes(self):
        raise NotImplementedError()

    def _world_pos(self):
        xyz = self.stat.mom1()[::-1]
        if self.wcs is not None:
            return self.wcs.all_pix2world([xyz], self.wcs_origin).ravel()[::-1]
        return xyz[::-1]

    def sky_maj(self):
        """Major axis of the projection onto the PP plane

        Intensity weighted second moment in direction
        of greatest elongation in the PP plane
        """
        dx = self.dx
        a, b = self._sky_paxes()
        return dx * np.sqrt(self.stat.mom2_along(a))

    def sky_min(self):
        """Minor axis of the projection onto the PP plane

        Intensity-weighted second moment perpendicular
        to major axis, in PP plane
        """
        dx = self.dx
        a, b = self._sky_paxes()
        return dx * np.sqrt(self.stat.mom2_along(b))

    def sky_radius(self):
        """ Geometric mean of sky_maj and sky_min """
        u, a = _qsplit(self.sky_maj())
        u, b = _qsplit(self.sky_min())
        return u * np.sqrt(a * b)

    def sky_deconvolved_rad(self):
        """sky_radius corrected for beam-smearing"""
        beam = self.bmaj * self.bmin
        u, a = _qsplit(self.sky_maj())
        u, b = _qsplit(self.sky_min())
        return u * np.sqrt(np.sqrt(a ** 2 - beam) * np.sqrt(b ** 2 - beam))


class PPVStatistic(SpatialBase):
    dv = MetaData('dv', 'Velocity channel width')
    vaxis = MetaData('vaxis', 'Index of velocity axis (numpy convention)')

    def __init__(self, stat, metadata):
        """
        Compute properties of structures in a PPV cube

        Parameters
        ----------
        stat :   ScalarStatistic instance
        metadata : dict
                 Key-value paris of metadata
        """
        self.stat = stat
        self.metadata = metadata

    def _sky_paxes(self):
        vaxis = self.vaxis
        ax = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        ax.pop(vaxis)
        a, b = self.stat.projected_paxes(ax)
        a = list(a)
        a.insert(0, vaxis)
        b = list(b)
        b.insert(0, vaxis)
        return a, b

    def xcen(self):
        p = self._world_pos()
        return p[2] if self.vaxis != 2 else p[1]

    def ycen(self):
        p = self._world_pos()
        return p[1] if self.vaxis == 0 else p[0]

    def vcen(self):
        p = self._world_pos()
        return p[self.vaxis]

    def flux(self):
        """Integrated flux

        sum(v_i * dx^2 * dv)
        """
        fac = self.bunit * self.dx ** 2 * self.dv
        return fac * self.stat.mom0()

    def vrms(self):
        """Intensity-weighted second moment of velocity"""
        ax = [0, 0, 0]
        ax[self.vaxis] = 1
        return self.dv * np.sqrt(self.stat.mom2_along(ax))

    def sky_pa(self):
        """The position angle of sky_maj, sky_min

        Returns the angle in degrees counter-clockwise from the +x axis
        """
        a, b = self._sky_paxes()
        a.pop(self.vaxis)
        return np.degrees(np.arctan2(a[0], a[1]))


class PPStatistic(SpatialBase):

    def __init__(self, stat, metadata):
        self.stat = stat
        self.metadata = metadata

    def _sky_paxes(self):
        return self.stat.paxes()

    def flux(self):
        """ Integrated flux """
        fac = self.bunit * self.dx ** 2
        return fac * self.stat.mom0()

    def sky_pa(self):
        """The position angle of sky_maj, sky_min

        Returns the angle in degrees counter-clockwise from the +x axis
        """
        a, b = self._sky_paxes()
        return np.degrees(np.arctan2(a[0], a[1]))

    def xcen(self):
        return self._world_pos()[1]

    def ycen(self):
        return self._world_pos()[0]



class PPPStatistic(object):

    def __init__(self, rhostat, vstat, metadata):
        """
        Derive properties from PPP density and velocity fields

        This is not currently implemented

        Parameters
        ----------
        rhostat : ScalarStatistic instance
        vstat : VectorStatistic instance
        """
        raise NotImplementedError()

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


def _make_catalog(structures, fields, metadata, statistic, verbose):
    _warn_missing_metadata(statistic, metadata, verbose=verbose)

    result = None

    for struct in structures:
        if isinstance(struct, Structure):
            stat = ScalarStatistic(struct.values(), struct.indices())
        else:
            stat = ScalarStatistic(struct.values, struct.indices)
        stat = statistic(stat, metadata)
        row = dict((lbl, getattr(stat, lbl)())
                   for lbl in fields)
        if hasattr(struct, 'idx'):
            row.update(_idx=struct.idx)

        # first row
        if result is None:
            result = Table(names=sorted(row.keys()))
            for k, v in row.items():
                result[k].units = _unit(v)

        result.add_row(row)

    return result


def ppv_catalog(structures, metadata, fields=None, verbose=True):
    """
    Iterate over a collection of position-position-velocity (PPV) structures,
    extracting several quantities from each, and building a catalog

    Parameters
    ----------
    structures : iterable of Structures
         The structures to catalog (e.g., a dendrogram)
    metadata : dict
        The metadata used to compute the catalog
    fields : list of strings, optional
        The quantities to extract. If not provided,
        defaults to all PPV statistics
    verbose : bool, optional
        If True (the default), will generate warnings
        about missing metadata

    Returns
    -------
    table : a :class:`~astropy.table.table.Table` instance
        The resulting catalog
    """
    fields = fields or ['flux', 'luminosity', 'sky_maj',
                        'sky_min', 'sky_radius', 'sky_deconvolved_rad',
                        'sky_pa', 'vrms', 'xcen', 'ycen', 'vcen']
    return _make_catalog(structures, fields, metadata, PPVStatistic, verbose)


def pp_catalog(structures, metadata, fields=None, verbose=False):
    """
    Iterate over a collection of position-position (PP) structures, extracting
    several quantities from each, and building a catalog.

    Parameters
    ----------
    structures : iterable of Structures
         The structures to catalog (e.g., a dendrogram)
    metadata : dict
        The metadata used to compute the catalog
    fields : list of strings, optional
        The quantities to extract. If not provided,
        defaults to all PPV statistics
    verbose : bool, optional
        If True (the default), will generate warnings
        about missing metadata

    Returns
    -------
    table : a :class:`~astropy.table.table.Table` instance
        The resulting catalog
    """
    fields = fields or ['flux', 'luminosity', 'sky_maj',
                        'sky_min', 'sky_radius', 'sky_deconvolved_rad',
                        'sky_pa', 'xcen', 'ycen']
    return _make_catalog(structures, fields, metadata, PPStatistic, verbose)

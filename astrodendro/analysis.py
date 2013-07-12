# Licensed under an MIT open source license - see LICENSE

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
    if q is None:
        return None
    elif isinstance(1 * q, Quantity):
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
            Location of each element of values. The i-th array in the tuple
            describes the ith positional dimension
        """
        self.values = values.astype(np.float)
        self.indices = indices

    def mom0(self):
        """The sum of the values"""
        return np.nansum(self.values)

    def mom1(self):
        """The intensity-weighted mean position"""
        m0 = self.mom0()
        return [np.nansum(i * self.values) / m0 for i in self.indices]

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
        warnings.warn("%s missing, defaulting to %s" %
                      (m, m.default))


class SpatialBase(object):

    spatial_scale = MetaData('spatial_scale', 'Angular length of a pixel')
    data_unit = MetaData('data_unit', 'Units of the pixel values')
    distance = MetaData('distance', 'Distance')
    wcs = MetaData('wcs', 'WCS object', default=None)

    def _sky_paxes(self):
        raise NotImplementedError()

    def _world_pos(self):
        xyz = self.stat.mom1()[::-1]
        if self.wcs is not None:
            # We use origin=0 since the indices come from Numpy indexing
            return self.wcs.all_pix2world([xyz], 0).ravel()[::-1]
        return xyz[::-1]

    @property
    def major_sigma(self):
        """
        Major axis of the projection onto the position-position (PP) plane,
        computed from the intensity weighted second moment in direction of
        greatest elongation in the PP plane.
        """
        dx = self.spatial_scale
        a, b = self._sky_paxes()
        # We need to multiply the second moment by two to get the major axis
        # rather than the half-major axis.
        return dx * np.sqrt(self.stat.mom2_along(a))

    @property
    def minor_sigma(self):
        """
        Minor axis of the projection onto the position-position (PP) plane,
        computed from the intensity weighted second moment perpendicular to
        the major axis in the PP plane.
        """
        dx = self.spatial_scale
        a, b = self._sky_paxes()
        # We need to multiply the second moment by two to get the minor axis
        # rather than the half-minor axis.
        return dx * np.sqrt(self.stat.mom2_along(b))

    @property
    def radius(self):
        """
        Geometric mean of major_sigma and minor_sigma.
        """
        u, a = _qsplit(self.major_sigma)
        u, b = _qsplit(self.minor_sigma)
        return u * np.sqrt(a * b)


class PPVStatistic(SpatialBase):
    """
    Compute properties of structures in a position-position-velocity (PPV)
    cube.

    Parameters
    ----------
    structure : `~astrodendro.structure.Structure` instance
        The structure to compute the statistics for
    metadata : dict
         Key-value pairs of metadata
    """

    velocity_scale = MetaData('velocity_scale', 'Velocity channel width')
    vaxis = MetaData('vaxis', 'Index of velocity axis (numpy convention)',
                     default=0)

    def __init__(self, stat, metadata=None):
        if isinstance(stat, Structure):
            self.stat = ScalarStatistic(stat.values(subtree=True),
                                        stat.indices(subtree=True))
        else:
            self.stat = stat
        self.metadata = metadata or {}

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

    @property
    def x_cen(self):
        """
        The mean position of the structure in the x direction.
        """
        p = self._world_pos()
        return p[2] if self.vaxis != 2 else p[1]

    @property
    def y_cen(self):
        """
        The mean position of the structure in the y direction.
        """
        p = self._world_pos()
        return p[1] if self.vaxis == 0 else p[0]

    @property
    def v_cen(self):
        """
        The mean velocity of the structure.
        """
        p = self._world_pos()
        return p[self.vaxis]

    @property
    def v_rms(self):
        """
        Intensity-weighted second moment of velocity
        """
        ax = [0, 0, 0]
        ax[self.vaxis] = 1
        return self.velocity_scale * np.sqrt(self.stat.mom2_along(ax))

    @property
    def position_angle(self):
        """
        The position angle of sky_maj, sky_min in degrees counter-clockwise
        from the +x axis.
        """
        a, b = self._sky_paxes()
        a.pop(self.vaxis)
        return np.degrees(np.arctan2(a[0], a[1]))


class PPStatistic(SpatialBase):
    """
    Compute properties of structures in a position-position (PP) cube.

    Parameters
    ----------
    structure : `~astrodendro.structure.Structure` instance
        The structure to compute the statistics for
    metadata : dict
         Key-value pairs of metadata
    """

    def __init__(self, stat, metadata=None):
        if isinstance(stat, Structure):
            self.stat = ScalarStatistic(stat.values(subtree=True),
                                        stat.indices(subtree=True))
        else:
            self.stat = stat
        self.metadata = metadata or {}

    def _sky_paxes(self):
        return self.stat.paxes()

    @property
    def position_angle(self):
        """
        The position angle of sky_maj, sky_min in degrees counter-clockwise
        from the +x axis.
        """
        a, b = self._sky_paxes()
        return np.degrees(np.arctan2(a[0], a[1]))

    @property
    def x_cen(self):
        """
        The mean position of the structure in the x direction.
        """
        return self._world_pos()[1]

    @property
    def y_cen(self):
        """
        The mean position of the structure in the y direction.
        """
        return self._world_pos()[0]


class PPPStatistic(object):

    def __init__(self, rhostat, vstat, metadata=None):
        """
        Derive properties from PPP density and velocity fields

        This is not currently implemented

        Parameters
        ----------
        rhostat : ScalarStatistic instance
        vstat : VectorStatistic instance
        """
        raise NotImplementedError()

    @property
    def mass(self):
        pass

    @property
    def volume(self):
        pass

    @property
    def surface_area(self):
        pass

    @property
    def virial(self):
        pass

    @property
    def v_rms(self):
        pass

    @property
    def vz_rms(self):
        pass

    @property
    def pressure_vz(self):
        pass

    @property
    def pressure(self):
        pass


def _make_catalog(structures, fields, metadata, statistic, verbose):
    """
    Make a catalog from a list of structures
    """

    _warn_missing_metadata(statistic, metadata, verbose=verbose)

    result = None

    for struct in structures:
        stat = ScalarStatistic(struct.values(subtree=True),
                               struct.indices(subtree=True))
        stat = statistic(stat, metadata)
        row = dict((lbl, getattr(stat, lbl))
                   for lbl in fields)
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
    fields = fields or ['major_sigma', 'minor_sigma', 'radius',
                        'position_angle', 'v_rms', 'x_cen', 'y_cen', 'v_cen']
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
    fields = fields or ['major_sigma', 'minor_sigma', 'radius',
                        'position_angle', 'x_cen', 'y_cen']
    return _make_catalog(structures, fields, metadata, PPStatistic, verbose)

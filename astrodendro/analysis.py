# Licensed under an MIT open source license - see LICENSE

import abc
import warnings
from functools import wraps
from weakref import WeakKeyDictionary

import numpy as np

from astropy.units import Quantity
from astropy.table import Table
from astropy import units as u
from astropy.wcs import WCS

from . import six
from .structure import Structure
from .flux import UnitMetadataWarning

__all__ = ['ppv_catalog', 'pp_catalog']


def memoize(func):

    # cache[instance][method args] -> method result
    # hold weakrefs to instances,
    # to stay out of the way of the garbage collector
    cache = WeakKeyDictionary()

    @wraps(func)
    def wrapper(self, *args):
        try:
            return cache[self][args]
        except KeyError:
            cache.setdefault(self, {})[args] = func(self, *args)
            return cache[self][args]
        except TypeError:
            warnings.warn("Cannot memoize inputs to %s" % func)
            return func(self, *args)

    return wrapper


class MissingMetadataWarning(UserWarning):
    pass


def _qsplit(q):
    """Split a potential astropy Quantity into unit/quantity"""
    if isinstance(1. * q, Quantity):
        return q.unit, q.value

    return 1, q


def _unit(q):
    """Return the units associated with a number, array, unit, or Quantity"""
    if q is None:
        return None
    elif isinstance(1 * q, Quantity):
        return (1 * q).unit


class ScalarStatistic(object):
    # This class does all of the heavy computation

    def __init__(self, values, indices):
        """
        Compute pixel-level statistics from a scalar field, sampled at specific
        locations.

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

    @memoize
    def mom0(self):
        """The sum of the values"""
        return np.nansum(self.values)

    @memoize
    def mom1(self):
        """The intensity-weighted mean position"""
        m0 = self.mom0()
        return [np.nansum(i * self.values) / m0 for i in self.indices]

    @memoize
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

    @memoize
    def mom2_along(self, direction):
        """
        Intensity-weighted variance/covariance along 1 or more directions.

        Parameters
        ----------
        direction : array like
                  One or more set of direction vectors. Need not be normalized

        Returns
        -------
        result : array
            The variance (or co-variance matrix) of the data along the
            specified direction(s).
        """
        w = np.atleast_2d(direction).astype(np.float)
        for row in w:
            row /= np.linalg.norm(row)

        result = np.dot(np.dot(w, self.mom2()), w.T)
        if result.size == 1:
            result = np.asscalar(result)
        return result

    @memoize
    def paxes(self):
        """
        The principal axes of the data (direction of greatest elongation)

        Returns
        -------
        result : tuple
            Ordered tuple of ndarrays

        Notes
        -----
        Each array is a normalized direction vector. The arrays
        are sorted in decreasing order of elongation of the data
        """
        mom2 = self.mom2()
        w, v = np.linalg.eig(mom2)
        order = np.argsort(w)

        return tuple(v[:, o] for o in order[::-1])

    @memoize
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
        result : tuple
            Tuple of arrays (nnew items)

        Notes
        -----
        The ordered principal axes in the new space
        """
        axes = tuple(axes)
        mom2 = self.mom2_along(axes)
        w, v = np.linalg.eig(mom2)
        order = np.argsort(w)

        return tuple(v[:, o] for o in order[::-1])

    @memoize
    def count(self):
        """
        Number of elements in the dataset.
        """
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


class Metadata(object):

    """
    A descriptor to wrap around metadata dictionaries.

    Lets classes reference self.x instead of self.metadata['x'],
    """

    _restrict_types = None

    def __init__(self, key, description, default=None, strict=False):
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
        if not isinstance(key, six.string_types):
            raise TypeError("Key is", key, type(key))
        self.key = key
        self.description = description or 'no description'
        self.default = default
        self.strict = strict

    def __get__(self, instance, type=None):

        if instance is None:
            return self

        try:
            value = instance.metadata[self.key]
        except KeyError:
            if self.strict:
                raise KeyError("Required metadata item not found: %s" % self)
            else:
                if self.default is not None:
                    warnings.warn("{0} ({1}) missing, defaulting to {2}".format(self.key, self.description, self.default),
                                  MissingMetadataWarning)
                value = self.default

        if value is not None and self._restrict_types is not None:
            if isinstance(value, self._restrict_types):
                return value
            else:
                raise TypeError("{0} should be an instance of {1}".format(self.key, ' or '.join([x.__name__ for x in self._restrict_types])))
        else:
            return value

    def __str__(self):
        return "%s (%s)" % (self.key, self.description)


class MetadataQuantity(Metadata):
    _restrict_types = (u.UnitBase, u.Quantity)


class MetadataWCS(Metadata):
    _restrict_types = (WCS,)


class SpatialBase(object):

    __metaclass__ = abc.ABCMeta

    wavelength = MetadataQuantity('wavelength', 'Wavelength')
    spatial_scale = MetadataQuantity('spatial_scale', 'Pixel width/height')
    beam_major = MetadataQuantity('beam_major', 'Major FWHM of beam')
    beam_minor = MetadataQuantity('beam_minor', 'Minor FWHM of beam')
    data_unit = MetadataQuantity('data_unit', 'Units of the pixel values', strict=True)
    wcs = MetadataWCS('wcs', 'WCS object')

    @abc.abstractmethod
    def _sky_paxes(self):
        raise NotImplementedError()

    def _world_pos(self):
        xyz = self.stat.mom1()[::-1]
        if self.wcs is None:
            return xyz[::-1] * u.pixel
        else:
            # TODO: set units correctly following WCS
            # We use origin=0 since the indices come from Numpy indexing
            return self.wcs.all_pix2world([xyz], 0).ravel()[::-1]

    @abc.abstractproperty
    def flux(self):
        raise NotImplementedError

    @abc.abstractproperty
    def x_cen(self):
        raise NotImplementedError()

    @abc.abstractproperty
    def y_cen(self):
        raise NotImplementedError()

    @abc.abstractproperty
    def position_angle(self):
        raise NotImplementedError()

    @property
    def major_sigma(self):
        """
        Major axis of the projection onto the position-position (PP) plane,
        computed from the intensity weighted second moment in direction of
        greatest elongation in the PP plane.
        """
        dx = self.spatial_scale or u.pixel
        a, b = self._sky_paxes()
        # We need to multiply the second moment by two to get the major axis
        # rather than the half-major axis.
        return dx * np.sqrt(self.stat.mom2_along(tuple(a)))

    @property
    def minor_sigma(self):
        """
        Minor axis of the projection onto the position-position (PP) plane,
        computed from the intensity weighted second moment perpendicular to
        the major axis in the PP plane.
        """
        dx = self.spatial_scale or u.pixel
        a, b = self._sky_paxes()
        # We need to multiply the second moment by two to get the minor axis
        # rather than the half-minor axis.
        return dx * np.sqrt(self.stat.mom2_along(tuple(b)))

    @property
    def radius(self):
        """
        Geometric mean of ``major_sigma`` and ``minor_sigma``.
        """
        u, a = _qsplit(self.major_sigma)
        u, b = _qsplit(self.minor_sigma)
        return u * np.sqrt(a * b)

    @property
    def area_ellipse(self):
        """
        The area of the ellipse defined by the second moments, where the
        semi-major and semi-minor axes used are the HWHM (half-width at
        half-maximum) derived from the moments.
        """
        return np.pi * self.major_sigma * self.minor_sigma * (2.3548 * 0.5) ** 2

    def to_mpl_ellipse(self, **kwargs):
        """
        Returns a Matplotlib ellipse representing the first and second moments
        of the structure.

        Any keyword arguments are passed to :class:`~matplotlib.patches.Ellipse`
        """
        from matplotlib.patches import Ellipse
        return Ellipse((self.x_cen.value, self.y_cen.value),
                       self.major_sigma.value * 2.3548,
                       self.minor_sigma.value * 2.3548,
                       angle=self.position_angle.value,
                       **kwargs)


class PPVStatistic(SpatialBase):

    """
    Compute properties of structures in a position-position-velocity (PPV)
    cube.

    Parameters
    ----------
    structure : :class:`~astrodendro.structure.Structure` instance
        The structure to compute the statistics for
    metadata : dict
         Key-value pairs of metadata
    """

    velocity_scale = MetadataQuantity('velocity_scale', 'Velocity channel width')
    vaxis = Metadata('vaxis', 'Index of velocity axis (numpy convention)', default=0)

    def __init__(self, stat, metadata=None):
        if isinstance(stat, Structure):
            self.stat = ScalarStatistic(stat.values(subtree=True),
                                        stat.indices(subtree=True))
        else:
            self.stat = stat
        if len(self.stat.indices) != 3:
            raise ValueError("PPVStatistic can only be used on 3-d datasets")
        self.metadata = metadata or {}

    def _sky_paxes(self):
        vaxis = self.vaxis
        ax = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        ax.pop(vaxis)
        a, b = self.stat.projected_paxes(tuple(ax))
        a = list(a)
        a.insert(0, vaxis)
        b = list(b)
        b.insert(0, vaxis)
        return tuple(a), tuple(b)

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
        The mean velocity of the structure (where the velocity axis can be
        specified by the ``vaxis`` metadata parameter, which defaults to 0
        following the Numpy convention - the third axis in the FITS convention).
        """
        p = self._world_pos()
        return p[self.vaxis]

    @property
    def flux(self):
        """
        The integrated flux of the structure, in Jy (note that this does not
        include any kind of background subtraction, and is just a plain sum of
        the values in the structure, converted to Jy).
        """
        from .flux import compute_flux
        return compute_flux(self.stat.mom0() * self.data_unit,
                            u.Jy,
                            wavelength=self.wavelength,
                            spatial_scale=self.spatial_scale,
                            velocity_scale=self.velocity_scale,
                            beam_major=self.beam_major,
                            beam_minor=self.beam_minor)

    @property
    def v_rms(self):
        """
        Intensity-weighted second moment of velocity (where the velocity axis
        can be specified by the ``vaxis`` metadata parameter, which defaults to
        0 following the Numpy convention - the third axis in the FITS
        convention).
        """
        dv = self.velocity_scale or u.pixel
        ax = [0, 0, 0]
        ax[self.vaxis] = 1
        return dv * np.sqrt(self.stat.mom2_along(tuple(ax)))

    @property
    def position_angle(self):
        """
        The position angle of sky_maj, sky_min in degrees counter-clockwise
        from the +x axis (note that this is the +x axis in pixel coordinates,
        which is the ``-x`` axis for conventional astronomy images).
        """
        a, b = self._sky_paxes()
        a = list(a)
        a.pop(self.vaxis)
        return np.degrees(np.arctan2(a[0], a[1])) * u.degree

    @property
    def area_exact(self):
        """
        The exact area of the structure on the sky.
        """
        dx = self.spatial_scale or u.pixel
        indices = zip(*tuple(self.stat.indices[i] for i in range(3) if i != self.vaxis))
        return len(set(indices)) * dx ** 2


class PPStatistic(SpatialBase):

    """
    Compute properties of structures in a position-position (PP) cube.

    Parameters
    ----------
    structure : :class:`~astrodendro.structure.Structure` instance
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
        if len(self.stat.indices) != 2:
            raise ValueError("PPStatistic can only be used on 2-d datasets")
        self.metadata = metadata or {}

    def _sky_paxes(self):
        return self.stat.paxes()

    @property
    def flux(self):
        """
        The integrated flux of the structure, in Jy (note that this does not
        include any kind of background subtraction, and is just a plain sum of
        the values in the structure, converted to Jy).
        """
        from .flux import compute_flux
        return compute_flux(self.stat.mom0() * self.data_unit,
                            u.Jy,
                            wavelength=self.wavelength,
                            spatial_scale=self.spatial_scale,
                            beam_major=self.beam_major,
                            beam_minor=self.beam_minor)

    @property
    def position_angle(self):
        """
        The position angle of sky_maj, sky_min in degrees counter-clockwise
        from the +x axis.
        """
        a, b = self._sky_paxes()
        return np.degrees(np.arctan2(a[0], a[1])) * u.degree

    @property
    def x_cen(self):
        """
        The mean position of the structure in the x direction (in pixel
        coordinates, or in world coordinates if the WCS transformation is
        available in the meta-data).
        """
        return self._world_pos()[1]

    @property
    def y_cen(self):
        """
        The mean position of the structure in the y direction (in pixel
        coordinates, or in world coordinates if the WCS transformation is
        available in the meta-data).
        """
        return self._world_pos()[0]

    @property
    def area_exact(self):
        """
        The exact area of the structure on the sky.
        """
        dx = self.spatial_scale or u.pixel
        return self.stat.count() * dx ** 2


class PPPStatistic(object):

    def __init__(self, rhostat, vstat, metadata=None):
        """
        Derive properties from PPP density and velocity fields.

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


def _make_catalog(structures, fields, metadata, statistic):
    """
    Make a catalog from a list of structures
    """

    result = None

    for struct in structures:
        stat = ScalarStatistic(struct.values(subtree=True),
                               struct.indices(subtree=True))
        stat = statistic(stat, metadata)
        row = {}
        for lbl in fields:
            row[lbl] = getattr(stat, lbl)

        row = dict((lbl, getattr(stat, lbl))
                   for lbl in fields)
        row.update(_idx=struct.idx)

        # first row
        if result is None:
            sorted_row_keys = sorted(row.keys())
            try:
                result = Table(names=sorted_row_keys,
                               dtype=[int if x == '_idx' else float for x in sorted_row_keys])
            except TypeError:  # dtype was called dtypes in older versions of Astropy
                result = Table(names=sorted_row_keys,
                               dtypes=[int if x == '_idx' else float for x in sorted_row_keys])
            for k, v in row.items():
                try:  # Astropy API change
                    result[k].unit = _unit(v)
                except AttributeError:
                    result[k].units = _unit(v)

        # astropy.table.Table should in future support setting row items from
        # quantities, but for now we need to strip off the quantities
        new_row = {}
        for x in row:
            if row[x] is not None:  # in Astropy 0.3+ we no longer need to exclude None items
                if isinstance(row[x], Quantity):
                    new_row[x] = row[x].value
                else:
                    new_row[x] = row[x]
        result.add_row(new_row)

    result.sort('_idx')

    return result


def ppv_catalog(structures, metadata, fields=None, verbose=True):
    """
    Iterate over a collection of position-position-velocity (PPV) structures,
    extracting several quantities from each, and building a catalog.

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
    fields = fields or ['major_sigma', 'minor_sigma', 'radius', 'area_ellipse', 'area_exact',
                        'position_angle', 'v_rms', 'x_cen', 'y_cen', 'v_cen', 'flux']
    with warnings.catch_warnings():
        warnings.simplefilter("once" if verbose else 'ignore', category=MissingMetadataWarning)
        warnings.simplefilter("once" if verbose else 'ignore', category=UnitMetadataWarning)        
        return _make_catalog(structures, fields, metadata, PPVStatistic)


def pp_catalog(structures, metadata, fields=None, verbose=True):
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
    fields = fields or ['major_sigma', 'minor_sigma', 'radius', 'area_ellipse', 'area_exact',
                        'position_angle', 'x_cen', 'y_cen', 'flux']
    with warnings.catch_warnings():
        warnings.simplefilter("once" if verbose else 'ignore', category=MissingMetadataWarning)
        return _make_catalog(structures, fields, metadata, PPStatistic)

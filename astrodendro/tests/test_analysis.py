# Licensed under an MIT open source license - see LICENSE

import pytest
from mock import patch

import numpy as np
from numpy.testing import assert_allclose
import astropy.units as u
from astropy.wcs import WCS


from ._testdata import data
from ..analysis import (ScalarStatistic, PPVStatistic, ppv_catalog,
                        Metadata, PPStatistic, pp_catalog)
from .. import Dendrogram
from ..structure import Structure


def assert_allclose_quantity(a, b):
    if not isinstance(a, u.Quantity):
        raise TypeError("a is not a quantity")
    if not isinstance(b, u.Quantity):
        raise TypeError("b is not a quantity")
    assert_allclose(a.value, b.value)
    assert a.unit == b.unit

wcs_2d = WCS(header=dict(cdelt1=1, crval1=0, crpix1=1,
                         cdelt2=2, crval2=0, crpix2=1))
wcs_3d = WCS(header=dict(cdelt1=1, crval1=0, crpix1=1,
                         cdelt2=2, crval2=0, crpix2=1,
                         cdelt3=3, crval3=0, crpix3=1))


def benchmark_stat():
    x = np.array([216, 216, 216, 216, 216, 217, 216,
                  216, 216, 217, 216, 217, 218, 216,
                  217, 216, 216, 217, 216, 217, 216])
    y = np.array([48, 50, 51, 47, 48, 48, 49, 50,
                  46, 46, 47, 47, 47, 48, 48, 49,
                  46, 46, 47, 46, 47])
    z = np.array([11, 11, 11, 12, 12, 12, 12, 12,
                  13, 13, 13, 13, 13, 13, 13, 13,
                  14, 14, 14, 15, 15])
    indices = (z, y, x)
    values = data[indices]
    return ScalarStatistic(values, indices)


def benchmark_values():
    result = {}
    result['mom0'] = 45.8145142
    result['mom1'] = [12.8093996, 47.6574211, 216.3799896]
    result['mom2'] = [
        [1.2831592097072804, -1.1429260442312184, 0.1653071908000249],
        [-1.1429260442312184, 2.0107706426889038, -0.2976682802749759],
        [0.1653071908000249, -0.2976682802749759, 0.3306051333187136]]
    result['mom2_100'] = 1.2831592097072804
    result['mom2_010'] = 2.0107706426889038
    result['mom2_011'] = 0.8730196376110217

    result['paxis1'] = [0.5833025351190873,
                        -0.8016446861607931, 0.1308585101314019]
    result['paxis2'] = [0.8038500408313988,
                        0.5466096691013617, -0.2346124069614796]
    result['paxis3'] = [-0.1165472624260408,
                        -0.2420406304632855, -0.9632409194100563]
    result['sig_maj'] = np.sqrt(2.0619485)
    result['sig_min'] = np.sqrt(0.27942730)
    result['area_exact_pp'] = 21
    result['area_exact_ppv'] = 10

    return result


class TestScalarStatistic(object):

    def setup_method(self, method):
        self.stat = benchmark_stat()

    def test_mom0(self):
        stat = self.stat
        assert_allclose(stat.mom0(), benchmark_values()['mom0'])

    def test_mom1(self):
        stat = self.stat
        assert_allclose(stat.mom1(), benchmark_values()['mom1'])

    def test_mom2(self):
        stat = self.stat
        assert_allclose(stat.mom2(), benchmark_values()['mom2'])

    def test_mom2_along(self):
        stat = self.stat
        v = benchmark_values()
        assert_allclose(stat.mom2_along((1, 0, 0)), v['mom2_100'])
        assert_allclose(stat.mom2_along((2, 0, 0)), v['mom2_100'])
        assert_allclose(stat.mom2_along((0, 1, 0)), v['mom2_010'])
        assert_allclose(stat.mom2_along((0, 1, 1)), v['mom2_011'])

    def test_count(self):
        stat = self.stat
        assert_allclose(stat.count(), 21)

    def test_sky_paxes(self):
        stat = self.stat
        v1, v2, v3 = stat.paxes()
        v = benchmark_values()

        # doesn't matter if v = vex * -1.
        assert_allclose(np.abs(np.dot(v1, v['paxis1'])), 1)
        assert_allclose(np.abs(np.dot(v2, v['paxis2'])), 1)
        assert_allclose(np.abs(np.dot(v3, v['paxis3'])), 1)

    def test_projected_paxes(self):
        stat = self.stat
        v = benchmark_values()
        v1, v2 = stat.projected_paxes(((0, 1, 0), (0, 0, 1)))

        assert_allclose(stat.mom2_along((0, v1[0], v1[1])), v['sig_maj'] ** 2)
        assert_allclose(stat.mom2_along((0, v2[0], v2[1])), v['sig_min'] ** 2)

    def test_projected_paxes_int(self):
        ind = (np.array([0, 1, 2]),
               np.array([0, 1, 2]),
               np.array([0, 1, 2]))
        v = np.array([1, 1, 1])
        stat = ScalarStatistic(v, ind)
        a, b = stat.projected_paxes(((0, 1, 0), (0, 0, 1)))
        assert_allclose(a, [1 / np.sqrt(2), 1 / np.sqrt(2)])


class TestScalar2D(object):

    def setup_method(self, method):
        x = np.array([213, 213, 214, 211, 212, 212])
        y = np.array([71, 71, 71, 71, 71, 71])
        z = np.array([46, 45, 45, 44, 46, 43])

        ind = (z, y, x)
        val = data[ind].copy()
        ind = (z, x)
        val[0] = 0  # we'll replace this with nan in a subclass

        self.stat = ScalarStatistic(val, ind)

    def test_mom0(self):
        assert_allclose(self.stat.mom0(), 19.2793083)

    def test_mom1(self):
        assert_allclose(self.stat.mom1(), [44.6037369, 212.4050598])

    def test_mom2(self):
        assert_allclose(self.stat.mom2(),
                        [[1.0321983103326871, 0.3604276691031663],
                        [0.3604276691031663, 1.0435691076146387]])

    def test_mom2_along(self):
        assert_allclose(self.stat.mom2_along((0, 1)), 1.0435691076146387)

    def test_count(self):
        assert_allclose(self.stat.count(), 6)

    def test_sky_paxes(self):
        v1, v2 = self.stat.paxes()
        v1ex = [0.7015083489012299, 0.7126612353859795]
        v2ex = [0.7126612353859795, -0.7015083489012299]

        assert_allclose(np.abs(np.dot(v1, v1ex)), 1)
        assert_allclose(np.abs(np.dot(v2, v2ex)), 1)


class TestScalarNan(TestScalar2D):

    def setup_method(self, method):
        x = np.array([213, 213, 214, 211, 212, 212])
        y = np.array([71, 71, 71, 71, 71, 71])
        z = np.array([46, 45, 45, 44, 46, 43])

        ind = (z, y, x)
        val = data[ind].copy()
        ind = (z, x)
        # all tests should be the same as superclass,
        # since nan should = 0 weight
        val[0] = np.nan

        self.stat = ScalarStatistic(val, ind)


class TestPPVStatistic(object):

    def setup_method(self, method):
        self.stat = benchmark_stat()
        self.v = benchmark_values()

    def metadata(self, **kwargs):
        result = dict(data_unit=u.Jy, wcs=wcs_3d)
        result.update(**kwargs)
        return result

    def test_x_cen(self):

        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.x_cen, self.v['mom1'][2])

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.x_cen, self.v['mom1'][1] * 2)

    def test_y_cen(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.y_cen, self.v['mom1'][1] * 2)

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.y_cen, self.v['mom1'][0] * 3)

    def test_v_cen(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.v_cen, self.v['mom1'][0] * 3)

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.v_cen, self.v['mom1'][2])

    def test_major_sigma(self):
        p = PPVStatistic(self.stat, self.metadata(spatial_scale=2 * u.arcsec))
        assert_allclose_quantity(p.major_sigma, self.v['sig_maj'] * 2 * u.arcsec)

    def test_minor_sigma(self):
        p = PPVStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.minor_sigma, self.v['sig_min'] * 4 * u.arcsec)

    def test_radius(self):
        p = PPVStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.radius, np.sqrt(self.v['sig_min'] * self.v['sig_maj']) * 4 * u.arcsec)

    def test_area_exact(self):
        p = PPVStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.area_exact, (4 * u.arcsec) ** 2 * self.v['area_exact_ppv'])

    def test_area_ellipse(self):
        p = PPVStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.area_ellipse, (4 * u.arcsec) ** 2 * self.v['sig_min'] * self.v['sig_maj'] * np.pi * (2.3548 * 0.5) ** 2)

    def test_v_rms(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose_quantity(p.v_rms, np.sqrt(self.v['mom2_100']) * u.pixel)

        p = PPVStatistic(self.stat, self.metadata(vaxis=1, velocity_scale=10 * u.km / u.s))
        assert_allclose_quantity(p.v_rms, np.sqrt(self.v['mom2_010']) * 10 * u.km / u.s)

    def test_position_angle(self):
        x = np.array([0, 1, 2])
        y = np.array([1, 1, 1])
        z = np.array([0, 1, 2])
        v = np.array([1, 1, 1])

        ind = (z, y, x)
        stat = ScalarStatistic(v, ind)
        p = PPVStatistic(stat, self.metadata())
        assert_allclose_quantity(p.position_angle, 0 * u.degree)

        ind = (z, x, y)
        stat = ScalarStatistic(v, ind)
        p = PPVStatistic(stat, self.metadata())
        assert_allclose_quantity(p.position_angle, 90 * u.degree)

    def test_units(self):
        m = self.metadata(spatial_scale=1 * u.deg, velocity_scale=1 * u.km / u.s,
                          data_unit=1 * u.K, distance=1 * u.pc)
        p = PPVStatistic(self.stat, m)

        assert p.v_rms.unit == u.km / u.s
        assert p.major_sigma.unit == u.deg
        assert p.minor_sigma.unit == u.deg
        assert p.radius.unit == u.deg


class TestPPStatistic(object):

    def setup_method(self, method):
        self.stat = benchmark_stat()
        # this trick essentially collapses along the 0th axis
        # should preserve sky_maj, sky_min
        self.stat.indices = (self.stat.indices[1], self.stat.indices[2])
        self.v = benchmark_values()

    def metadata(self, **kwargs):
        result = dict()
        result.update(**kwargs)
        return result

    def test_major_sigma(self):
        p = PPStatistic(self.stat, self.metadata(spatial_scale=2 * u.arcsec))

        assert_allclose_quantity(p.major_sigma, self.v['sig_maj'] * 2 * u.arcsec)

    def test_minor_sigma(self):
        p = PPStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.minor_sigma, self.v['sig_min'] * 4 * u.arcsec)

    def test_radius(self):
        p = PPStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.radius, np.sqrt(self.v['sig_min'] * self.v['sig_maj']) * 4 * u.arcsec)

    def test_area_exact(self):
        p = PPStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.area_exact, (4 * u.arcsec) ** 2 * self.v['area_exact_pp'])

    def test_area_ellipse(self):
        p = PPStatistic(self.stat, self.metadata(spatial_scale=4 * u.arcsec))
        assert_allclose_quantity(p.area_ellipse, (4 * u.arcsec) ** 2 * self.v['sig_min'] * self.v['sig_maj'] * np.pi * (2.3548 * 0.5) ** 2)

    def test_position_angle(self):
        x = np.array([0, 1, 2])
        y = np.array([1, 1, 1])
        v = np.array([1, 1, 1])

        ind = (y, x)
        stat = ScalarStatistic(v, ind)
        p = PPStatistic(stat, self.metadata())
        assert_allclose_quantity(p.position_angle, 0 * u.degree)

        ind = (x, y)
        stat = ScalarStatistic(v, ind)
        p = PPStatistic(stat, self.metadata())
        assert_allclose_quantity(p.position_angle, 90 * u.degree)


def test_statistic_dimensionality():

    d = Dendrogram.compute(np.ones((10, 10)))

    with pytest.raises(ValueError) as exc:
        PPVStatistic(d.trunk[0])
    assert exc.value.args[0] == "PPVStatistic can only be used on 3-d datasets"

    PPStatistic(d.trunk[0])

    d = Dendrogram.compute(np.ones((10, 10, 10)))

    with pytest.raises(ValueError) as exc:
        PPStatistic(d.trunk[0])
    assert exc.value.args[0] == "PPStatistic can only be used on 2-d datasets"

    PPVStatistic(d.trunk[0])


class TestCataloger(object):
    files = []
    cataloger = None

    def stat(self):
        raise NotImplementedError()

    def metadata(self):
        raise NotImplementedError()

    def make_catalog(self, s=None, md=None, fields=None):
        s = s or [self.stat()]
        structures = [Structure(zip(*x.indices), x.values) for x in s]
        md = md or self.metadata()
        return self.cataloger(structures, md, fields)

    def test_benchmark(self):
        c = self.make_catalog()
        assert c.dtype.names == tuple(sorted(self.fields))
        assert len(c) == 1
        c = self.make_catalog(s=[self.stat()] * 3)
        assert len(c) == 3
        return c

    def test_field_selection(self):
        stat = self.stat()
        md = self.metadata()
        c = self.cataloger([Structure(zip(*stat.indices), stat.values)], md, fields=['x_cen'])
        assert c.dtype.names == ('_idx', 'x_cen',)


class TestPPVCataloger(TestCataloger):
    fields = ['_idx', 'flux',
              'major_sigma', 'minor_sigma', 'radius', 'area_ellipse',
              'area_exact', 'v_rms', 'position_angle', 'x_cen', 'y_cen', 'v_cen']
    cataloger = staticmethod(ppv_catalog)

    def stat(self):
        return benchmark_stat()

    def metadata(self):
        return dict(vaxis=1, data_unit=u.Jy, wcs=wcs_3d)


class TestPPCataloger(TestCataloger):
    fields = ['_idx', 'flux',
              'major_sigma', 'minor_sigma', 'radius', 'area_ellipse',
              'area_exact', 'position_angle', 'x_cen', 'y_cen']
    cataloger = staticmethod(pp_catalog)

    def stat(self):
        bs = benchmark_stat()
        bs.indices = (bs.indices[1], bs.indices[2])
        return bs

    def metadata(self):
        return dict(data_unit=u.Jy, wcs=wcs_2d)


# don't let pytest test abstract class
del TestCataloger


def test_metadata_protocol():
    class Foo(object):
        x = Metadata('x', 'test')
        y = Metadata('y', 'test', default=5)
        z = Metadata('z', 'test', strict=True)

        def __init__(self, md):
            self.metadata = md

    f = Foo(dict(x=10))
    assert f.x == 10
    assert f.y == 5
    with pytest.raises(KeyError):
        f.z

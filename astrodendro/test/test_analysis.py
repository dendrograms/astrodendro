import pytest
from mock import patch

import numpy as np
from numpy.testing import assert_allclose
import astropy.units as u
from astropy.wcs import WCS


from ._testdata import data
from ..analysis import (ScalarStatistic, PPVStatistic, ppv_catalog,
                        _missing_metadata, MetaData, _warn_missing_metadata,
                        PPStatistic, pp_catalog)
from .. import Dendrogram
from ..structure import Structure


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
        assert_allclose(stat.mom2_along([1, 0, 0]), v['mom2_100'])
        assert_allclose(stat.mom2_along([2, 0, 0]), v['mom2_100'])
        assert_allclose(stat.mom2_along([0, 1, 0]), v['mom2_010'])
        assert_allclose(stat.mom2_along([0, 1, 1]), v['mom2_011'])

    def test_count(self):
        stat = self.stat
        assert_allclose(stat.count(), 21)

    def test_paxes(self):
        stat = self.stat
        v1, v2, v3 = stat.paxes()
        v = benchmark_values()

        #doesn't matter if v = vex * -1.
        assert_allclose(np.abs(np.dot(v1, v['paxis1'])), 1)
        assert_allclose(np.abs(np.dot(v2, v['paxis2'])), 1)
        assert_allclose(np.abs(np.dot(v3, v['paxis3'])), 1)

    def test_projected_paxes(self):
        stat = self.stat
        v = benchmark_values()
        v1, v2 = stat.projected_paxes([[0, 1, 0], [0, 0, 1]])

        assert_allclose(stat.mom2_along([0, v1[0], v1[1]]), v['sig_maj'] ** 2)
        assert_allclose(stat.mom2_along([0, v2[0], v2[1]]), v['sig_min'] ** 2)

    def test_projected_paxes_int(self):
        ind = (np.array([0, 1, 2]),
               np.array([0, 1, 2]),
               np.array([0, 1, 2]))
        v = np.array([1, 1, 1])
        stat = ScalarStatistic(v, ind)
        a, b = stat.projected_paxes([[0, 1, 0], [0, 0, 1]])
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
        assert_allclose(self.stat.mom2_along([0, 1]), 1.0435691076146387)

    def test_count(self):
        assert_allclose(self.stat.count(), 6)

    def test_paxes(self):
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
        #all tests should be the same as superclass,
        #since nan should = 0 weight
        val[0] = np.nan

        self.stat = ScalarStatistic(val, ind)


class TestPPVStatistic(object):
    def setup_method(self, method):
        self.stat = benchmark_stat()
        self.v = benchmark_values()

    def metadata(self, **kwargs):
        result = dict(distance=1,
                      dx=1,
                      dv=1,
                      vaxis=0,
                      bunit=1,
                      wcs=wcs_3d,
                      )
        result.update(**kwargs)
        return result

    def test_flux(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.flux(), self.v['mom0'])

        p = PPVStatistic(self.stat, self.metadata(dx=5))
        assert_allclose(p.flux(), self.v['mom0'] * 25)

        p = PPVStatistic(self.stat, self.metadata(dv=3))
        assert_allclose(p.flux(), self.v['mom0'] * 3)

    def test_xcen(self):

        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.xcen(), self.v['mom1'][2])

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.xcen(), self.v['mom1'][1] * 2)

    def test_ycen(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.ycen(), self.v['mom1'][1] * 2)

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.ycen(), self.v['mom1'][0] * 3)

    def test_vcen(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.vcen(), self.v['mom1'][0] * 3)

        p = PPVStatistic(self.stat, self.metadata(vaxis=2))
        assert_allclose(p.vcen(), self.v['mom1'][2])

    def test_sky_major_sigma(self):
        p = PPVStatistic(self.stat, self.metadata(dx=2))

        assert_allclose(p.sky_major_sigma(), self.v['sig_maj'] * 2)

    def test_sky_minor_sigma(self):
        p = PPVStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_minor_sigma(), self.v['sig_min'] * 4)

    def test_sky_radius(self):
        p = PPVStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_radius(), np.sqrt(self.v['sig_min'] *
                                                self.v['sig_maj']) * 4)

    def test_sky_vrms(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.vrms(), np.sqrt(self.v['mom2_100']))

        p = PPVStatistic(self.stat, self.metadata(vaxis=1, dv=10))
        assert_allclose(p.vrms(), np.sqrt(self.v['mom2_010']) * 10)

    def test_pa(self):
        x = np.array([0, 1, 2])
        y = np.array([1, 1, 1])
        z = np.array([0, 1, 2])
        v = np.array([1, 1, 1])

        ind = (z, y, x)
        stat = ScalarStatistic(v, ind)
        p = PPVStatistic(stat, self.metadata())
        assert_allclose(p.sky_pa(), 0)

        ind = (z, x, y)
        stat = ScalarStatistic(v, ind)
        p = PPVStatistic(stat, self.metadata())
        assert_allclose(p.sky_pa(), 90)

    def test_deconvolved_rad(self):
        p = PPVStatistic(self.stat, self.metadata(bmaj=.4, bmin=.1))

        a = self.v['sig_maj']
        b = self.v['sig_min']
        dcr = np.sqrt(np.sqrt(a ** 2 - .04) * np.sqrt(b ** 2 - .04))
        assert_allclose(p.sky_deconvolved_rad(), dcr)

    def test_luminosity(self):
        p = PPVStatistic(self.stat, self.metadata(dist=10))
        v = benchmark_values()
        assert_allclose(p.luminosity(), v['mom0'] * 100 * np.radians(1) ** 2)

        p = PPVStatistic(self.stat, self.metadata(dist=10, dx=1 * u.rad))
        assert_allclose(p.luminosity(), v['mom0'] * 100)

    def test_units(self):
        m = self.metadata(dx=1 * u.deg, dv=1 * u.km / u.s,
                          bunit=1 * u.K, dist=1 * u.pc)
        p = PPVStatistic(self.stat, m)

        assert p.vrms().unit == u.km / u.s
        assert p.flux().unit == u.deg ** 2 * u.km / u.s * u.K
        assert p.sky_major_sigma().unit == u.deg
        assert p.sky_minor_sigma().unit == u.deg
        assert p.sky_radius().unit == u.deg
        assert p.luminosity().unit == u.km / u.s * u.K * u.pc ** 2


class TestPPStatistic(object):
    def setup_method(self, method):
        self.stat = benchmark_stat()
        #this trick essentially collapses along the 0th axis
        #should preserve sky_maj, sky_min
        self.stat.indices = (self.stat.indices[1], self.stat.indices[2])
        self.v = benchmark_values()

    def metadata(self, **kwargs):
        result = dict(distance=1,
                      dx=1,
                      )
        result.update(**kwargs)
        return result

    def test_flux(self):
        p = PPStatistic(self.stat, self.metadata(dx=5))
        assert_allclose(p.flux(), self.v['mom0'] * 25)

    def test_sky_major_sigma(self):
        p = PPStatistic(self.stat, self.metadata(dx=2))

        assert_allclose(p.sky_major_sigma(), self.v['sig_maj'] * 2)

    def test_sky_minor_sigma(self):
        p = PPStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_minor_sigma(), self.v['sig_min'] * 4)

    def test_sky_radius(self):
        p = PPStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_radius(), np.sqrt(self.v['sig_min'] *
                                                self.v['sig_maj']) * 4)

    def test_pa(self):
        x = np.array([0, 1, 2])
        y = np.array([1, 1, 1])
        v = np.array([1, 1, 1])

        ind = (y, x)
        stat = ScalarStatistic(v, ind)
        p = PPStatistic(stat, self.metadata())
        assert_allclose(p.sky_pa(), 0)

        ind = (x, y)
        stat = ScalarStatistic(v, ind)
        p = PPStatistic(stat, self.metadata())
        assert_allclose(p.sky_pa(), 90)


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
        c = ppv_catalog([Structure(zip(*stat.indices), stat.values)], md, fields=['flux'])
        assert c.dtype.names == ('_idx', 'flux',)


class TestPPVCataloger(TestCataloger):
    fields = ['_idx', 'flux', 'luminosity',
              'sky_major_sigma', 'sky_minor_sigma', 'sky_radius',
              'vrms', 'sky_deconvolved_rad',
              'sky_pa', 'xcen', 'ycen', 'vcen']
    cataloger = staticmethod(ppv_catalog)

    def stat(self):
        return benchmark_stat()

    def metadata(self):
        return dict(vaxis=1, dx=1, dv=1, dist=1, lum2mass=1,
                    bmaj=1, bmin=1, bunit=1, wcs=wcs_3d)


class TestPPCataloger(TestCataloger):
    fields = ['_idx', 'flux', 'luminosity',
              'sky_major_sigma', 'sky_minor_sigma', 'sky_radius',
              'sky_deconvolved_rad',
              'sky_pa', 'xcen', 'ycen']
    cataloger = staticmethod(pp_catalog)

    def stat(self):
        bs = benchmark_stat()
        bs.indices = (bs.indices[1], bs.indices[2])
        return bs

    def metadata(self):
        return dict(dx=1, dist=1, lum2mass=1,
                    bmaj=1, bmin=1, bunit=1, wcs=wcs_2d)


#don't let pytest test abstract class
del TestCataloger


def test_find_missing_ppv_metadata():
    md = dict(dx=1, dv=1, vaxis=1, bmaj=1, bmin=1, bunit=1, dist=1,
              wcs=wcs_3d)
    assert len(_missing_metadata(PPVStatistic, md)) == 0

    md.pop('dx')
    assert _missing_metadata(PPVStatistic, md)[0].key == 'dx'
    assert len(_missing_metadata(PPVStatistic, {})) == 8


def test_metadata_protocol():
    class Foo(object):
        x = MetaData('x', 'test')
        y = MetaData('y', 'test', default=5)
        z = MetaData('z', 'test', strict=True)

        def __init__(self, md):
            self.metadata = md

    f = Foo(dict(x=10))
    assert f.x == 10
    assert f.y == 5
    with pytest.raises(KeyError):
        f.z


def test_warn_missing_metadata():
    class Foo(object):
        x = MetaData('x', 'test description')

    class Bar(object):
        y = MetaData('y', 'test', strict=True)

    with patch('warnings.warn') as mock:
        _warn_missing_metadata(Foo, {'x': 3})
    assert mock.call_count == 0

    with patch('warnings.warn') as mock:
        _warn_missing_metadata(Foo, {})
    assert mock.call_count == 1

    with patch('warnings.warn') as mock:
        _warn_missing_metadata(Foo, {}, verbose=False)
    assert mock.call_count == 0

    with pytest.raises(RuntimeError):
        _warn_missing_metadata(Bar, {})


def test_dendrogram_ppv_catalog():
    x = np.random.random((5, 5, 5))
    d = Dendrogram.compute(x)
    c = ppv_catalog(d, {})
    for ct, st in zip(c['flux'], d):
        assert ct == st.values(subtree=True).sum()


def test_dendrogram_ppv_catalog():
    x = np.random.random((5, 5))
    d = Dendrogram.compute(x)
    c = pp_catalog(d, {})
    for ct, st in zip(c['flux'], d):
        assert ct == st.values(subtree=True).sum()

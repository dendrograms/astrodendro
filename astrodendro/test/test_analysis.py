import numpy as np
from numpy.testing import assert_allclose
import astropy.units as u

from ._testdata import data
from ..analysis import ScalarStatistic, PPVStatistic, ppv_catalog

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
                        [0.3604276691031663,  1.0435691076146387]])

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
                      )
        result.update(**kwargs)
        return result

    def test_flux(self):
        p = PPVStatistic(self.stat, self.metadata())
        assert_allclose(p.flux(),  self.v['mom0'])

        p = PPVStatistic(self.stat, self.metadata(dx=5))
        assert_allclose(p.flux(),  self.v['mom0'] * 25)

        p = PPVStatistic(self.stat, self.metadata(dv=3))
        assert_allclose(p.flux(), self.v['mom0'] * 3)

    def test_sky_maj(self):
        p = PPVStatistic(self.stat, self.metadata(dx=2))

        assert_allclose(p.sky_maj(), self.v['sig_maj'] * 2)

    def test_sky_min(self):
        p = PPVStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_min(), self.v['sig_min'] * 4)

    def test_sky_rad(self):
        p = PPVStatistic(self.stat, self.metadata(dx=4))
        assert_allclose(p.sky_rad(), np.sqrt(self.v['sig_min'] *
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
        assert_allclose(p.luminosity(), v['mom0'] * 100 * np.radians(1)**2)

        p = PPVStatistic(self.stat, self.metadata(dist=10, dx=1 * u.rad))
        assert_allclose(p.luminosity(), v['mom0'] * 100)

    def test_units(self):
        m = self.metadata(dx = 1 * u.deg, dv = 1 * u.km / u.s,
                          bunit = 1 * u.K, dist = 1 * u.pc)
        p = PPVStatistic(self.stat, m)

        assert p.vrms().unit == u.km / u.s
        assert p.flux().unit == u.deg ** 2 * u.km / u.s * u.K
        assert p.sky_maj().unit == u.deg
        assert p.sky_min().unit == u.deg
        assert p.sky_rad().unit == u.deg
        assert p.luminosity().unit == u.km / u.s * u.K * u.pc ** 2


class TestPPVCataloger(object):
    def test_benchmark(self):
        stat = benchmark_stat()
        #anything with .values and .indices will do
        s = [stat]

        md = dict(vaxis=1, dx=1, dv=1, dist=1, lum2mass=1,
                  bmaj=1, bmin=1, bunit=1)
        c = ppv_catalog(s, md)
        assert len(c) == 1
        assert set(c[0].keys()) == set(['flux', 'luminosity',
                                        'sky_maj', 'sky_min', 'sky_rad',
                                        'vrms', 'sky_deconvolved_rad',
                                        'sky_pa'])

        s = [stat, stat, stat]
        c = ppv_catalog(s, md)
        assert len(c) == 3


    def test_field_selection(self):
        stat = benchmark_stat()
        md = dict(vaxis=1, dx=1, dv=1, dist=1, lum2mass=1,
                  bmaj=1, bmin=1, bunit=1)
        c = ppv_catalog([stat], md, fields=['flux'])
        assert c[0].keys() == ['flux']

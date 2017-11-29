import pytest

import fixture
import route
from math import *


class PythonImpl(object):
    EARTH_RADIUS = 6356752.0
    WGS84_a = 6378137.0
    WGS84_b = 6356752.3142
    WGS84_f = 1 / 298.257223563

    @staticmethod
    def hav(x):
        return (sin(x)/2)**2

    @staticmethod
    def to_rad(*xs):
        return [x*pi/180 for x in xs]

    @classmethod
    def rlat(cls, x):
        return atan((1-cls.WGS84_f) * tan(x))

    @classmethod
    def harversine(cls, lat1, lng1, lat2, lng2):
        lat1, lng1, lat2, lng2 = cls.to_rad(lat1, lng1, lat2, lng2)
        return 2 * cls.EARTH_RADIUS * asin(
            sqrt(cls.hav(lat2-lat1) + cos(lat1)*cos(lat2)*cls.hav(lng2-lng1)))

    @classmethod
    def vincenty(cls, lat1, lng1, lat2, lng2, limit=20):
        if lng1 == lng2 and lat1 == lat2:
            return .0
        lat1, lng1, lat2, lng2 = cls.to_rad(lat1, lng1, lat2, lng2)
        L = lng2 - lng1
        U1, U2 = cls.rlat(lat1), cls.rlat(lat2)
        sinU1, cosU1 = sin(U1), cos(U1)
        sinU2, cosU2 = sin(U2), cos(U2)
        sigma, sinSigma, cosSigma = 0, 0, 0
        sinAlpha, cosSqAlpha, cos2SigmaM = 0, 0, 0
        C, Lambda, LambdaP = 0, L, 2*pi

        while abs(Lambda - LambdaP) > 1e-12 and limit > 0:
            sinLambda, cosLambda = sin(Lambda), cos(Lambda)
            sinSigma = sqrt((cosU2*sinLambda)**2 +
                            (cosU1*sinU2-sinU1*cosU2*cosLambda)**2)
            cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
            sigma = atan2(sinSigma, cosSigma)
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
            cosSqAlpha = 1 - sinAlpha**2
            cos2SigmaM = (.0 if abs(cosSqAlpha) < 1e-12 else
                          cosSigma - 2*sinU1*sinU2/cosSqAlpha)
            C = cls.WGS84_f/16 * cosSqAlpha * (
                4 + cls.WGS84_f * (4 - 3 * cosSqAlpha))
            LambdaP = Lambda
            Lambda = L + (1-C) * cls.WGS84_f * sinAlpha * (
                sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM**2)))
            limit -= 1

        if limit == 0:
            return float('nan')

        uSq = cosSqAlpha * (cls.WGS84_a**2 - cls.WGS84_b**2) / cls.WGS84_b**2
        A = 1 + uSq/16384 * (4096 + uSq * (-768 + uSq * (320-175*uSq)))
        B = uSq/1024 * (256 + uSq * (-128 + uSq * (74-47*uSq)))
        deltaSigma = B * sinSigma * (cos2SigmaM +
                                     B/4 * cosSigma*(-1+2*cos2SigmaM**2) -
                                     B/6 * cos2SigmaM*(-3+4*sinSigma**2) *
                                     (-3+4*cos2SigmaM**2))

        return cls.WGS84_b * A * (sigma - deltaSigma)

    @classmethod
    def measure(cls, ps, formula):
        ds = [.0] * len(ps)
        for i in range(1, len(ps)):
            ds[i] = ds[i-1] + formula(ps[i-1][0], ps[i-1][1],
                                      ps[i][0], ps[i][1])
        return ds


def test_py_vincenty():
    assert abs(PythonImpl.vincenty(
        43.65331, -79.38277, 43.65244, -79.38243) - 100.478) < 0.1
    assert abs(PythonImpl.vincenty(0.001, 0, 0, 0) - 110.574) < 0.01
    assert abs(PythonImpl.vincenty(0, 0.001, 0, 0) - 111.319) < 0.01
    assert PythonImpl.vincenty(0, 0, 0, 0) == 0


def test_py_haversine():
    assert abs(PythonImpl.harversine(
        43.65331, -79.38277, 43.65244, -79.38243) - 100.478) < 1.0
    assert abs(PythonImpl.harversine(0.001, 0, 0, 0) - 110.574) < 1.0
    assert (PythonImpl.harversine(0.001, 0, 0, 0)
            ==
            PythonImpl.harversine(0, 0.001, 0, 0))
    assert PythonImpl.harversine(0, 0, 0, 0) == 0


@pytest.fixture(scope="session")
def points():
    return fixture.routes[0]['coordinates']


@pytest.mark.benchmark(group="distance")
def test_distance_vincenty_c(benchmark):
    benchmark(route.distance_vincenty,
              43.65331, -79.38277,
              43.65244, -79.38243)


@pytest.mark.benchmark(group="distance")
def test_distance_vincenty_py(benchmark):
    benchmark(PythonImpl.vincenty,
              43.65331, -79.38277,
              43.65244, -79.38243)


@pytest.mark.benchmark(group="distance")
def test_distance_haversine_c(benchmark):
    benchmark(route.distance_haversine,
              43.65331, -79.38277,
              43.65244, -79.38243)


@pytest.mark.benchmark(group="distance")
def test_distance_haversine_py(benchmark):
    benchmark(PythonImpl.harversine,
              43.65331, -79.38277,
              43.65244, -79.38243)


@pytest.mark.benchmark(group="measure")
def test_measure_vincenty_c(benchmark, points):
    benchmark(route.measure_vincenty, points)


@pytest.mark.benchmark(group="measure")
def test_measure_vincenty_py(benchmark, points):
    benchmark(PythonImpl.measure, points, PythonImpl.vincenty)


@pytest.mark.benchmark(group="measure")
def test_measure_haversine_c(benchmark, points):
    benchmark(route.measure_haversine, points)


@pytest.mark.benchmark(group="measure")
def test_measure_haversine_py(benchmark, points):
    benchmark(PythonImpl.measure, points, PythonImpl.harversine)

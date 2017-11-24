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

    @classmethod
    def rlat(cls, x):
        return atan((1-cls.WGS84_f) * tan(x))

    @classmethod
    def harversine(cls, lng1, lat1, lng2, lat2):
        return 2 * cls.EARTH_RADIUS * asin(
            sqrt(cls.hav(lat2-lat1) + cos(lat1)*cos(lat2)*cls.hav(lng2-lng1)))

    @classmethod
    def vincenty(cls, lng1, lat1, lng2, lat2, limit=20):
        if lng1 == lng2 and lat1 == lat2:
            return .0
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


@pytest.fixture(scope="session")
def points():
    return map(lambda t: (t[0]*pi/180, t[1]*pi/180),
               fixture.routes[0]['coordinates'])


@pytest.mark.benchmark(group="distance")
def test_distance_vincenty_c(benchmark):
    benchmark(route.distance,
              -1.385490705853379,
              0.7618939888937657,
              -1.291657616475959,
              0.7105718841071718)


@pytest.mark.benchmark(group="distance")
def test_distance_vincenty_py(benchmark):
    benchmark(PythonImpl.vincenty,
              -1.385490705853379,
              0.7618939888937657,
              -1.291657616475959,
              0.7105718841071718)


@pytest.mark.benchmark(group="distance")
def test_distance_haversine_c(benchmark):
    benchmark(route.distance,
              -1.385490705853379,
              0.7618939888937657,
              -1.291657616475959,
              0.7105718841071718,
              route.HAVERSINE)


@pytest.mark.benchmark(group="distance")
def test_distance_haversine_py(benchmark):
    benchmark(PythonImpl.harversine,
              -1.385490705853379,
              0.7618939888937657,
              -1.291657616475959,
              0.7105718841071718)


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

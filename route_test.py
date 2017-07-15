import math
import route
import fixture


def test_measure_vincenty():
    for r in fixture.routes:
        ps = [(p[0]*math.pi/180, p[1]*math.pi/180) for p in r['coordinates']]
        ds1 = route.measure(ps, route.VINCENTY)
        assert abs(ds1[-1] - r['distance']) < 5, r['name']


def test_measure_haversine():
    for r in fixture.routes:
        ps = [(p[0]*math.pi/180, p[1]*math.pi/180) for p in r['coordinates']]
        ds1 = route.measure(ps, route.HAVERSINE)
        assert abs(ds1[-1] - r['distance']) < 5000, r['name']


def test_distance_vincenty():
    assert abs(route.distance(1e-6, 0, 0, 0) - 6.378) < 0.01
    assert abs(route.distance(0, 1e-6, 0, 0) - 6.335) < 0.01
    assert route.distance(0, 0, 0, 0) == 0


def test_distance_haversine():
    assert abs(route.distance(1e-6, 0, 0, 0, route.HAVERSINE) - 6.357) < 0.01
    assert (route.distance(1e-6, 0, 0, 0, route.HAVERSINE)
            ==
            route.distance(0, 1e-6, 0, 0, route.HAVERSINE))

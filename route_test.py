import route
import fixture


def test_measure_vincenty():
    for r in fixture.routes:
        ds = route.measure_vincenty(r['coordinates'])
        assert abs(ds[-1] - r['distance']) < 5, r['name']


def test_measure_haversine():
    for r in fixture.routes:
        ds = route.measure_haversine(r['coordinates'])
        assert abs(ds[-1] - r['distance']) < 5000, r['name']


def test_distance_vincenty():
    assert abs(route.distance_vincenty(
        43.65331, -79.38277, 43.65244, -79.38243) - 100.478) < 0.1
    assert abs(route.distance_vincenty(0.001, 0, 0, 0) - 110.574) < 0.01
    assert abs(route.distance_vincenty(0, 0.001, 0, 0) - 111.319) < 0.01
    assert route.distance_vincenty(0, 0, 0, 0) == 0


def test_distance_haversine():
    assert abs(route.distance_haversine(
        43.65331, -79.38277, 43.65244, -79.38243) - 100.478) < 1.0
    assert abs(route.distance_haversine(0.001, 0, 0, 0) - 110.574) < 1.0
    assert (route.distance_haversine(0.001, 0, 0, 0)
            ==
            route.distance_haversine(0, 0.001, 0, 0))
    assert route.distance_haversine(0, 0, 0, 0) == 0

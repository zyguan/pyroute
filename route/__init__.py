import _route

VINCENTY = 1
HAVERSINE = 2


def distance(lat1, lng1, lat2, lng2, formula=VINCENTY, iterlimit=20):
    """
    Calculate distance between two points.
    :param lat1: latitude of point 1 in degrees
    :param lng1: longitude of point 1 in degrees
    :param lat2: latitude of point 2 in degrees
    :param lng2: longitude of point 2 in degrees
    :param formula: distance formula to use
    :param iterlimit: max iterations, used by vincenty's formula
    :return: the distance between (lat1, lng1) and (lat2, lng2) in meter
    """
    return _route.distance(lat1, lng1, lat2, lng2, formula, iterlimit)


def distance_vincenty(lat1, lng1, lat2, lng2):
    return _route.distance(lat1, lng1, lat2, lng2, VINCENTY)


def distance_haversine(lat1, lng1, lat2, lng2):
    return _route.distance(lat1, lng1, lat2, lng2, HAVERSINE)


def measure(ps, formula=VINCENTY, iterlimit=20):
    """
    Calculate cumulative distances of a route.
    :param ps: points in route, eg: [(lng1, lat1), (lng2, lat2), ...]
    :param formula: distance formula to use
    :param iterlimit: max iterations, used by vincenty's formula
    :return: the cumulative distance array `ds` where ds[j]-ds[i] is the
    distance from (lat_i, lng_i) to (lat_j, lng_j) in meter.
    """
    return _route.measure(ps, formula, iterlimit)


def measure_vincenty(ps):
    return _route.measure(ps, VINCENTY)


def measure_haversine(ps):
    return _route.measure(ps, HAVERSINE)

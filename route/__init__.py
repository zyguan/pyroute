import _route

VINCENTY = 1
HAVERSINE = 2


def distance(lng1, lat1, lng2, lat2, formula=VINCENTY, iterlimit=20):
    """
    Calculate distance between two points.
    :param lng1: longitude of point 1 in radians
    :param lat1: latitude of point 1 in radians
    :param lng2: longitude of point 2 in radians
    :param lat2: latitude of point 2 in radians
    :param formula: distance formula to use
    :param iterlimit: max iterations, used by vincenty's formula
    :return: the distance between (lng1, lat1) and (lng2, lat2)
    """
    return _route.distance(lng1, lat1, lng2, lat2, formula, iterlimit)


def measure(ps, formula=VINCENTY, iterlimit=20):
    """
    Calculate cumulative distances of a route.
    :param ps: points in route, eg: [(lng1, lat1), (lng2, lat2), ...]
    :param formula: distance formula to use
    :param iterlimit: max iterations, used by vincenty's formula
    :return: the cumulative distance array `ds` where ds[j]-ds[i] is the
    distance from (lng_i, lat_i) to (lng_j, lat_j).
    """
    return _route.measure(ps, formula, iterlimit)

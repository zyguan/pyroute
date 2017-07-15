#include <Python.h>

#define EARTH_RADIUS 6356752.0
#define WGS84_a 6378137.0
#define WGS84_b 6356752.3142
#define WGS84_f (1/298.257223563)
#define _square(x) ((x)*(x))
#define _hav(x) _square(sin((x)/2.0))

#define F_VINCENTY  1
#define F_HAVERSINE 2

static double
vincenty_distance(double lng1, double lat1, double lng2, double lat2, long iterLimit)
{
    if (lng1 == lng2 && lat1 == lat2) {
        return 0;
    }
    double L = lng2 - lng1;
    double U1 = atan((1-WGS84_f) * tan(lat1));
    double U2 = atan((1-WGS84_f) * tan(lat2));
    double sinU1 = sin(U1), cosU1 = cos(U1);
    double sinU2 = sin(U2), cosU2 = cos(U2);

    double sinSigma = 0;
    double cosSigma = 0;
    double sigma = 0;
    double sinAlpha = 0;
    double cosSqAlpha = 0;
    double cos2SigmaM = 0;
    double C = 0;
    double lambda = L, lambdaP = 2*M_PI;

    while (fabs(lambda - lambdaP) > 1e-12 && iterLimit > 0) {
        double sinLambda = sin(lambda);
        double cosLambda = cos(lambda);
        sinSigma = sqrt(_square(cosU2*sinLambda) + _square(cosU1*sinU2-sinU1*cosU2*cosLambda));
        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        sinAlpha = cosU1*cosU2*sinLambda / sinSigma;
        cosSqAlpha = 1 - _square(sinAlpha);
        if (fabs(cosSqAlpha) < 1e-12) {
            cos2SigmaM = 0;
        } else {
            cos2SigmaM = cosSigma - 2*sinU1*sinU2 / cosSqAlpha;
        }
        C = WGS84_f/16 * cosSqAlpha * (4 + WGS84_f*(4-3*cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1-C) * WGS84_f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C*cosSigma*(-1+2*_square(cos2SigmaM))));
        iterLimit--;
    }

    if (iterLimit == 0) {
        return NAN;
    }

    double uSq = cosSqAlpha * (_square(WGS84_a) - _square(WGS84_b)) / _square(WGS84_b);
    double A = 1 + uSq/16384 * (4096 + uSq * (-768 + uSq * (320-175*uSq)));
    double B = uSq/1024 * (256 + uSq * (-128 + uSq * (74-47*uSq)));
    double deltaSigma = B * sinSigma * (cos2SigmaM + B/4 * (cosSigma*(-1 + 2*_square(cos2SigmaM))) - B/6 * cos2SigmaM * (-3 + 4*_square(sinSigma)) * (-3 + 4*_square(cos2SigmaM)));

    return WGS84_b * A * (sigma-deltaSigma);
}

static double
haversine_distance(double lng1, double lat1, double lng2, double lat2)
{
    return 2.0 * EARTH_RADIUS * asin(sqrt(_hav(lat2-lat1) + cos(lat1)*cos(lat2)*_hav(lng2-lng1)));
}

static PyObject *
distance(PyObject *self, PyObject *args) {
    double lng1, lat1, lng2, lat2, d = 0;
    long limit = -1;
    long formula = -1;

    if (!PyArg_ParseTuple(args, "dddd|ll:distance",
                          &lng1, &lat1, &lng2, &lat2, &formula, &limit)) {
        return NULL;
    };

    if (formula < 0) {
        formula = F_VINCENTY;
    }
    if (formula != F_VINCENTY && formula != F_HAVERSINE) {
        PyErr_SetString(PyExc_ValueError, "invalid formula code");
        return NULL;
    }

    if (limit < 0) {
        limit = 20;
    }

    if (formula == F_VINCENTY) {
        d = vincenty_distance(lng1, lat1, lng2, lat2, limit);
    } else if (formula == F_HAVERSINE) {
        d = haversine_distance(lng1, lat1, lng2, lat2);
    }

    return PyFloat_FromDouble(d);
}

static PyObject *
measure(PyObject *self, PyObject *args) {
    PyObject *ps;
    long limit = -1;
    long formula = -1;

    if (!PyArg_ParseTuple(args, "O|ll:measure", &ps, &formula, &limit)) {
        return NULL;
    }

    if (formula < 0) {
        formula = F_VINCENTY;
    }
    if (formula != F_VINCENTY && formula != F_HAVERSINE) {
        PyErr_SetString(PyExc_ValueError, "invalid formula code");
        return NULL;
    }

    if (limit < 0) {
        limit = 20;
    }

    // check route type
    if (!PySequence_Check(ps)) {
        PyErr_SetString(PyExc_TypeError, "route must be a sequence of tuple(lng:float, lat:float)");
        return NULL;
    }

    // calc distances
    Py_ssize_t len = PySequence_Size(ps);
    PyObject *ds = PyList_New(len);
    if (ds == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate a new list for distances");
        return NULL;
    }

    double lng1, lat1, lng2, lat2, dsum = 0;

    for (Py_ssize_t i = 0; i < len; i++) {
        // check type
        PyObject *p = PySequence_GetItem(ps, i);
        if (!PyTuple_Check(p)) {
            Py_XDECREF(p); Py_DECREF(ds);
            PyErr_SetString(PyExc_TypeError, "route must be a sequence of tuple(lng:float, lat:float)");
            return NULL;
        }
        if (PyTuple_Size(p) != 2) {
            Py_XDECREF(p); Py_DECREF(ds);
            PyErr_SetString(PyExc_TypeError, "route must be a sequence of tuple(lng:float, lat:float)");
            return NULL;
        }
        PyObject *lng = PyTuple_GetItem(p, 0);
        if (!PyFloat_Check(lng)) {
            Py_XDECREF(p); Py_DECREF(ds);
            PyErr_SetString(PyExc_TypeError, "route must be a sequence of tuple(lng:float, lat:float)");
            return NULL;
        }
        PyObject *lat = PyTuple_GetItem(p, 1);
        if (!PyFloat_Check(lat)) {
            Py_XDECREF(p); Py_DECREF(ds);
            PyErr_SetString(PyExc_TypeError, "route must be a sequence of tuple(lng:float, lat:float)");
            return NULL;
        }

        lng2 = PyFloat_AsDouble(lng);
        lat2 = PyFloat_AsDouble(lat);
        Py_DECREF(p);

        if (i > 0) {
            if (formula == F_VINCENTY) {
                dsum += vincenty_distance(lng1, lat1, lng2, lat2, limit);
            } else if (formula == F_HAVERSINE) {
                dsum += haversine_distance(lng1, lat1, lng2, lat2);
            }
        }

        PyList_SetItem(ds, i, PyFloat_FromDouble(dsum));

        lng1 = lng2;
        lat1 = lat2;
    }

    return ds;
}

static PyMethodDef methods[] = {
    {"distance", distance, METH_VARARGS, "calculate distance of two geo points."},
    {"measure", measure, METH_VARARGS, "calculate cumulative distances of a route."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_route",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC
PyInit__route(void) {
    return PyModule_Create(&moduledef);
}

#else

PyMODINIT_FUNC
init_route(void) {
  (void) Py_InitModule("_route", methods);
}

#endif

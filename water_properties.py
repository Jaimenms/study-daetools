from iapws.iapws97 import _Region1, _TSat_P
from iapws._iapws import _ThCond, _Viscosity
import numpy as np
import inspect

def simple_model_parameters():

    c_dict = {
        'conductivity': np.array([[-6.70232734e-01, 3.72674263e-09, 1.71778072e-14],
                            [6.97259473e-03, -1.91154051e-11, -1.06811308e-16],
                            [-9.02411509e-06, 2.85075484e-14, 1.65519545e-19]]),
        'density': np.array([[7.78456318e+02, 2.40153974e-06, 8.29162678e-13],
                       [1.74126580e+00, -1.20391370e-08, -5.01742419e-15],
                       [-3.38166607e-03, 1.84233265e-11, 7.58071678e-18]]),
        'enthalpy': np.array([[-1.13185854e+06, 1.05561614e-03, 1.17467908e-08],
                        [4.11870305e+03, 1.20950026e-06, -7.29010857e-11],
                        [9.79273111e-02, -5.55494252e-09, 1.12754536e-13]]),
        'heat_capacity': np.array([[5.34437955e+03, -4.56546130e-05, -1.71806886e-12],
                             [-7.40001244e+00, 2.54359095e-07, 1.07068599e-14],
                             [1.17465273e-02, -3.72077253e-10, -1.64404045e-17]]),
        'viscosity': np.array([[1.89901972e-02, -2.70460526e-11, -4.75198595e-18],
                         [-1.03966136e-04, 1.60710675e-13, 2.94337924e-20],
                         [1.45088884e-07, -2.36355370e-16, -4.54320112e-23]]),
    }

    c_dict = {
 'conductivity': np.array([[-6.68466623e-01,  5.60495502e-09, -9.47189764e-17],
       [ 6.96267417e-03, -3.08101787e-11,  5.75613760e-19],
       [-9.01035342e-06,  4.66554489e-14, -8.77655266e-22]]),
 'density': np.array([[ 7.79090743e+02,  2.48303614e-06, -1.72867793e-14],
       [ 1.73778137e+00, -1.25464465e-08,  1.06919245e-16],
       [-3.37696769e-03,  1.92123247e-11, -1.69246313e-19]]),
 'enthalpy': np.array([[-1.13179941e+06,  2.33827218e-03, -2.06860245e-11],
       [ 4.11846511e+03, -6.75337539e-06,  1.19007291e-13],
       [ 9.80972572e-02,  6.76532230e-09, -1.71417931e-16]]),
 'heat_capacity': np.array([[ 5.33854120e+03, -4.59622909e-05,  6.92912051e-13],
       [-7.36540443e+00,  2.56423946e-07, -4.16832999e-15],
       [ 1.16954320e-02, -3.75497945e-10,  6.31920830e-18]]),
 'viscosity': np.array([[ 1.88609725e-02, -2.70404394e-11,  1.05272925e-19],
       [-1.03228185e-04,  1.60840574e-13, -6.17146995e-22],
       [ 1.44044684e-07, -2.36804252e-16,  9.04800488e-25]])}

    return c_dict


def simplified_model(x, y):

    c_dict = simple_model_parameters()

    c = c_dict[inspect.stack()[1][3]]

    f = 0.

    for i in range(c.shape[0]):

        for j in range(c.shape[1]):
            f += c[i][j] * (x ** i) * (y ** j)

    return f


def saturation_temperature(P):

    P = P*1e-6

    T = _TSat_P(P)

    if T > 623.15:

        return None

    else:

        return T


def density(T, P, simplified = False):

    if simplified:

        return simplified_model(T, P)

    P = P*1e-6

    v = _Region1(T, P)['v']
    alfav = _Region1(T, P)['alfav']
    kt = _Region1(T, P)['kt']

    rho = 1 / v # kg/m3
    drhodT = -alfav / v  # kg/m3K
    drhodP = kt / v * 1e-6 # kg/m3Pa

    return rho, drhodT, drhodP


def enthalpy(T, P, simplified = False):

    if simplified:

        return simplified_model(T, P)

    P = P*1e-6
    h = 1e+3 * _Region1(T, P)['h']
    cp = 1e+3 * _Region1(T, P)['cp']
    v = _Region1(T, P)['v']
    alfav = _Region1(T, P)['alfav']

    dvdT = alfav * v

    dhdP = v - T * dvdT
    dhdT = cp

    return h, dhdT, dhdP


def heat_capacity(T, P, simplified = False):

    if simplified:

        return simplified_model(T, P)

    P = P*1e-6
    cp = 1e+3 * _Region1(T, P)['cp']

    return cp, 0, 0


def conductivity(T, P, simplified = False):

    if simplified:

        return simplified_model(T, P)

    P = P*1e-6

    v = _Region1(T, P)['v']

    rho = 1 / v

    return _ThCond(rho, T), 0, 0


def viscosity(T, P, simplified = False):

    if simplified:

        return simplified_model(T, P)

    P = P*1e-6

    v = _Region1(T, P)['v']

    rho = 1 / v

    mu = _Viscosity(rho, T)

    return _Viscosity(rho, T), 0, 0


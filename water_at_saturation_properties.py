from iapws.iapws97 import _Region4, _Region1
from iapws._iapws import _ThCond, _Viscosity
import numpy as np
import inspect

def simple_model_parameters():

    c_dict = {
         'density': np.array([ 9.96906358e+02, -8.11519963e-04,  7.28615211e-09, -3.02080951e-14]),
         'conductivity': np.array([6.20940155e-01, 1.90022194e-06, -2.45476728e-11, 1.12611502e-16]),
         'enthalpy': np.array([1.31488947e+05, 7.28807358e+00, -7.84449119e-05, 3.44521411e-10]),
         'heat_capacity': np.array([4.17429701e+03, 4.31727962e-04, 7.44863129e-10, -8.92950974e-15]),
         'viscosity': np.array([6.97318785e-04, -1.47128076e-08, 1.98579652e-13, -9.38312502e-19])}

    return c_dict

def simplified_model(x,):

    c_dict = simple_model_parameters()

    c = c_dict[inspect.stack()[1][3]]

    f = 0.

    for i in range(c.shape[0]):
        f += c[i] * (x ** i)

    return f


def calculate_saturation_temperature(P):

    P = P*1e-6
    out = _Region4(P, 0.)

    return out['T']


def density(P, simplified = False):

    if simplified:

        return simplified_model(P)


    T = calculate_saturation_temperature(P)

    P = P*1e-6
    out = _Region1(T, P)

    v = out['v']

    rho = 1 / v # kg/m3

    return rho


def enthalpy(P, simplified = False):

    if simplified:

        return simplified_model(P)

    T = calculate_saturation_temperature(P)

    P = P*1e-6
    out = _Region1(T, P)
    h = 1e3 * out['h']

    return h


def heat_capacity(P, simplified = False):

    if simplified:

        return simplified_model(P)

    T = calculate_saturation_temperature(P)

    P = P*1e-6
    out = _Region1(T, P)
    cp = 1e3 * out['cp']

    return cp


def conductivity(P, simplified = False):

    if simplified:

        return simplified_model(P)

    T = calculate_saturation_temperature(P)

    P = P*1e-6
    out = _Region1(T, P)
    T = out['T']
    v = out['v']

    rho = 1 / v

    return _ThCond(rho, T)


def viscosity(P, simplified = False):

    if simplified:

        return simplified_model(P)

    T = calculate_saturation_temperature(P)

    P = P*1e-6
    out = _Region1(T, P)
    T = out['T']
    v = out['v']

    rho = 1 / v

    return _Viscosity(rho, T)

__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from pyUnits import m, kg, s, K, Pa, J, W, rad


class FixedExternalConvection(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        daeModelExtended.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)

    def define_parameters(self):

        self.Do = daeParameter("D_o", m, self, "Outside pipe diameter")
        self.kwall = daeParameter("k_wall", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-1)), self, "Wall conductivity")
        self.Text = daeParameter("T_ext", K, self, "External Temperature")
        self.hext = daeParameter("h_ext", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-2)), self, "External heat transfer coefficient")


    def define_variables(self):

        # Variable types
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10, 0.1, 1e-05)

        # Secondary variables
        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", [self.x, ])

        # State variables
        self.Tw = daeVariable("Tw", water_temperature_t, self, "Wall Temperature", [self.x, ])


    def eq_heat_balance(self):

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)

        A = 0.25 * 3.14 * self.D(x) ** 2
        eq.Residual = self.rho(x) * self.cp(x) * dt(A * self.T(x)) + self.k() * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)
        print("HEAT BALANCE FROM FIXED CONVECTION")


    def eq_total_he(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        prandtl = self.cp(x) * self.mu(x) / self.kappa(x)
        # nusselt =  0.027 * (Re ** (4/5)) * (prandtl ** (1/3))
        nusselt = (self.fD(x) / 8.) * (self.Re(x) - 1000.) * prandtl / (
                    1. + 12.7 * Sqrt(self.fD(x) / 8.) * (prandtl ** 2 / 3) - 1.)
        hint = nusselt * self.kappa(x) / self.D(x)
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.Di() * hext)
        Resint = 1 / (2 * self.pi * self.D(x) * hint)
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / (Resint + Reswall + Resext)


    def eq_wall_he(self):

        eq = self.CreateEquation("WallHeat", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.Do() * hext)
        eq.Residual = self.Qout(x) - (self.Tw(x) - self.Text()) / Resext


    def DeclareEquations(self):
        daeModelExtended.DeclareEquations(self)

        self.eq_heat_balance()
        self.eq_total_he()
        self.eq_wall_he()

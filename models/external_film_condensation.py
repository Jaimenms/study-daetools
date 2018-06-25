__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from pyUnits import m, kg, s, K, Pa, J, W, rad

from water_properties import density, viscosity, conductivity, heat_capacity

#from steam_properties import density, viscosity, conductivity, heat_capacity


class ExternalFilmCondensation(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        daeModelExtended.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)

    def define_parameters(self):

        self.Do = daeParameter("Do", m, self, "Outside pipe diameter")
        self.kwall = daeParameter("kwall", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-1)), self, "Wall conductivity")
        self.rhov = daeParameter("rhov", (kg / m**3), self, "External steam density")
        self.hvap = daeParameter("hvap", J / kg, self, "Vaporization heat")
        self.Tsat = daeParameter("Tsat", K, self, "Saturation temperature")
        self.Psat = daeParameter("Psat", Pa, self, "Saturation pressure")
        self.fNtub = daeParameter("fNtub", unit(), self, "Factor for number of pipes over the actual pipe, including it")


    def define_variables(self):

        # Variable types
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10, 0.1, 1e-05)
        thermal_resistance_t = daeVariableType("thermal_resistance_t", (K ** (1))*(W ** (-1))*(m ** (1)), -1e-10, 1e+10, 100, 1e-05)
        # heat_transfer_coefficient_t = daeVariableType("heat_transfer_coefficient_t",
        #                                               (K ** (-1)) * (W ** (1)) * (m ** (-2)), 0,
        #                                               1e+10, 10000, 1e-01)

        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", [self.x, ])

        self.To = daeVariable("To", water_temperature_t, self, "Outside Wall Temperature", [self.x, ])

        self.Ti = daeVariable("Ti", water_temperature_t, self, "Inside Wall Temperature", [self.x, ])

        self.hint = daeVariable("hint", heat_transfer_coefficient_t, self, "Internal Convection coefficient", [self.x, ])

        self.hext = daeVariable("hext", heat_transfer_coefficient_t, self, "External Convection coefficient", [self.x, ])

        self.Resistance = daeVariable("Resistance", thermal_resistance_t, self, "Overall Thermal Resistance", [self.x, ])


    def eq_heat_balance(self):

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)

        A = 0.25 * 3.14 * self.D(x) ** 2
        eq.Residual = self.rho(x) * self.cp(x) * dt(A * self.T(x)) + self.k() * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)


    def eq_calculate_hint(self):

        eq = self.CreateEquation("Hint", "Heat balance - Hint")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)

        # Calculates the Nussel dimensionless number using Petukhov correlation modified by Gnielinski. See Incropera 4th Edition [8.63]
        prandtl = self.cp(x) * self.mu(x) / self.kappa(x)
        nusselt = (self.fD(x) / 8.) * (self.Re(x) - 1000.) * prandtl / (
                1. + 12.7 * Sqrt(self.fD(x) / 8.) * (prandtl ** 2 / 3) - 1.)
        hint = nusselt * self.kappa(x) / self.D(x)

        eq.Residual = self.hint(x) - hint


    def eq_calculate_hext(self):

        eq = self.CreateEquation("Hext", "Heat balance - Hext")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)

        Tf = 0.5 * (self.Tsat() + self.To(x))
        Pf = self.Psat()

        rho_o = density(Tf / Constant(1 * K), Pf / Constant(1 * Pa), simplified=True) * Constant(1 * (kg ** (1))*(m ** (-3)))
        kappa_o = conductivity(Tf / Constant(1 * K), Pf / Constant(1 * Pa), simplified=True) * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))
        mu_o = viscosity(Tf / Constant(1 * K), Pf / Constant(1 * Pa), simplified=True) * Constant(1 * (Pa ** (1))*(s ** (1)))

        num = (self.g * rho_o * (rho_o - self.rhov()) * kappa_o ** 3 * self.hvap())
        den = mu_o * Abs(self.Tsat() - self.To(x)) * self.Do()

        hd1 = 0.729 * (num / den) ** 0.25

        eq.Residual = self.hext(x) - self.fNtub() * hd1


    def eq_calculate_resistance(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Resext = 1 / (2 * self.pi * self.Do() * self.hext(x))
        Resint = 1 / (2 * self.pi * self.D(x) * self.hint(x))
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        # TODO - Lembrar de colocar o Refilme no caso com Biofilme
        #Resfilm = Log(self.Di() / self.D()) / (2 * self.pi * self.kappa())
        eq.Residual = self.Resistance(x) - (Resint + Reswall + Resext)


    def eq_total_he(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) * self.Resistance(x) - (self.T(x) - self.Tsat())


    def eq_calculate_To(self):

        eq = self.CreateEquation("WallHeat1", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) / (2 * self.pi * self.Do() * self.hext(x)) - (self.To(x) - self.Tsat())


    def eq_calculate_Ti(self):

        eq = self.CreateEquation("WallHeat2", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) *  Log(self.Do() / self.Di()) - (self.Ti(x) - self.To(x)) * (2 * self.pi * self.kwall())


    def DeclareEquations(self):
        daeModelExtended.DeclareEquations(self)

        self.eq_heat_balance()
        self.eq_total_he()
        self.eq_wall_he()
        self.eq_calculate_resistance()

__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from pyUnits import m, kg, s, K, Pa, J, W, rad


class FixedExternalTemperature(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        daeModelExtended.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)

    def define_parameters(self):

        self.Do = daeParameter("Do", m, self, "Outside pipe diameter")
        self.kwall = daeParameter("kwall", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-1)), self, "Wall conductivity")
        self.Text = daeParameter("Text", K, self, "External Temperature")


    def define_variables(self):

        # Variable types
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10, 0.1, 1e-05)
        thermal_resistance_t = daeVariableType("thermal_resistance_t", (K ** (1))*(W ** (-1))*(m ** (1)), -1e-10, 1e+10, 100, 1e-05)

        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", [self.x, ])

        self.To = daeVariable("To", water_temperature_t, self, "Outside Wall Temperature", [self.x, ])

        self.Ti = daeVariable("Ti", water_temperature_t, self, "Internal Wall Temperature", [self.x, ])

        self.hint = daeVariable("hint", heat_transfer_coefficient_t, self, "Internal convection coefficient", [self.x, ])

        self.Resistance = daeVariable("Resistance", thermal_resistance_t, self, "Overall Thermal Resistance", [self.x, ])


    def eq_heat_balance(self):

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)

        A = 0.25 * 3.14 * self.D(x) ** 2
        eq.Residual = self.rho(x) * self.cp(x) * dt(A * self.T(x)) + self.k() * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)


    def eq_calculate_hint(self):

        eq = self.CreateEquation("InternalConvection", "Internal convection - hint")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)

        # Calculates the Nussel dimensionless number using Petukhov correlation modified by Gnielinski. See Incropera 4th Edition [8.63]
        prandtl = self.cp(x) * self.mu(x) / self.kappa(x)
        nusselt = (self.fD(x) / 8.) * (self.Re(x) - 1000.) * prandtl / (
                1. + 12.7 * Sqrt(self.fD(x) / 8.) * (prandtl ** 2 / 3) - 1.)
        hint = nusselt * self.kappa(x) / self.D(x)

        eq.Residual = self.hint(x) - hint
        #eq.Residual = self.hint(x) - 1398.8091902015963 * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-2)))


    def eq_calculate_resistance(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Resint = 1 / (2 * self.pi * self.D(x) * self.hint(x))
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        # TODO - Lembrar de colocar o Refilme no caso com Biofilme
        #Resfilm = Log(self.Di() / self.D()) / (2 * self.pi * self.kappa())
        eq.Residual = self.Resistance(x) - (Resint + Reswall)


    def eq_total_he(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x)*self.Resistance(x) - (self.T(x) - self.Text())


    def eq_calculate_To(self):

        eq = self.CreateEquation("WallHeat", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.To(x) - self.Text()


    def eq_calculate_Ti(self):

        eq = self.CreateEquation("WallHeat", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        eq.Residual = self.Qout(x) - (self.Ti(x) - self.To(x)) / Reswall


    def DeclareEquations(self):
        daeModelExtended.DeclareEquations(self)

        self.eq_heat_balance()
        self.eq_total_he()
        self.eq_calculate_To()
        self.eq_calculate_Ti()
        self.eq_calculate_hint()
        self.eq_calculate_resistance()
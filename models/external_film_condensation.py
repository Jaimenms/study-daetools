__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from pyUnits import m, kg, s, K, Pa, J, W, rad

from water_properties import density, viscosity, conductivity, heat_capacity
from daetools_extended.tools import daeVariable_wrapper, distribute_on_domains


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
        self.fNtub = daeParameter("fNtub", unit(), self, "Factor for number of pipes over the actual pipe, including it", self.YDomains)


    def define_variables(self):

        # Variable types
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10, 0.1, 1e-05)
        thermal_resistance_t = daeVariableType("thermal_resistance_t", (K ** (1))*(W ** (-1))*(m ** (1)), 1e-6, 100., 1e-6, 1e-05)
        heat_transfer_coefficient_t = daeVariableType("heat_transfer_coefficient_t",
                                                      (K ** (-1)) * (W ** (1)) * (m ** (-2)), 0.01,
                                                      1000000, 10000, 1e-01)

        self.To = daeVariable("To", water_temperature_t, self, "Outside Wall Temperature", self.Domains)
        self.Ti = daeVariable("Ti", water_temperature_t, self, "Inside Wall Temperature", self.Domains)
        self.hint = daeVariable("hint", heat_transfer_coefficient_t, self, "Internal Convection coefficient", self.Domains)
        self.hext = daeVariable("hext", heat_transfer_coefficient_t, self, "External Convection coefficient", self.Domains)
        self.Resistance = daeVariable("Resistance", thermal_resistance_t, self, "Overall Thermal Resistance", self.Domains)


    def eq_calculate_hint(self):

        eq = self.CreateEquation("InternalConvection", "Internal convection - hint")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        T = daeVariable_wrapper(self.T, domains)
        P = daeVariable_wrapper(self.P, domains)
        Ti = daeVariable_wrapper(self.Ti, domains)
        fD = daeVariable_wrapper(self.fD, domains)
        Re = daeVariable_wrapper(self.Re, domains)
        hint = daeVariable_wrapper(self.hint, domains)
        D = daeVariable_wrapper(self.D, domains)

        # Calculates the Nussel dimensionless number using Petukhov correlation modified by Gnielinski. See Incropera 4th Edition [8.63]
        Tm = 0.5 * (T + Ti)
        Tast = Tm / Constant(1 * K)
        Past = P / Constant(1 * Pa)

        mu = viscosity( Tast, Past, simplified = True)
        kappa = conductivity( Tast, Past, simplified = True)
        cp = heat_capacity( Tast, Past, simplified = True)
        prandtl = cp * mu / kappa
        kappa_i = kappa * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))

        nusselt = (fD / 8.) * (Re - 1000.) * prandtl / (1. + 12.7 * Sqrt(fD / 8.) * (prandtl ** (2 / 3)) - 1.)
        hint_calc = nusselt * kappa_i / D

        eq.Residual = hint - hint_calc


    def eq_calculate_hext(self):

        eq = self.CreateEquation("Hext", "Heat balance - Hext")
        xdomains = distribute_on_domains(self.XDomains, eq, eClosedClosed)
        ydomains = distribute_on_domains(self.YDomains, eq, eClosedClosed)
        domains = xdomains + ydomains

        g = self.g
        Tsat = self.Tsat()
        Psat = self.Psat()
        rhov = self.rhov()
        hvap = self.hvap()
        Do = self.Do()

        To = daeVariable_wrapper(self.To, domains)
        hext = daeVariable_wrapper(self.hext, domains)
        fNtub = daeVariable_wrapper(self.fNtub, ydomains)

        Tf = 0.5 * (Tsat + To)
        Pf = Psat

        Tast = Tf / Constant(1 * K)
        Past = Pf / Constant(1 * Pa)

        rho_o = density( Tast, Past, simplified = True) * Constant(1 * (kg ** (1))*(m ** (-3)))
        mu_o = viscosity( Tast, Past, simplified = True) * Constant(1 * (Pa ** (1))*(s ** (1)))
        kappa_o = conductivity( Tast, Past, simplified = True)  * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))

        num = (g * rho_o * (rho_o - rhov) * kappa_o ** 3 * hvap)

        den = mu_o * Abs(Tsat - To) * Do

        hd1 = (0.729 * (num / den) ** 0.25)

        eq.Residual = hext - fNtub * hd1


    def eq_calculate_resistance(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        hext = daeVariable_wrapper(self.hext, domains)
        hint = daeVariable_wrapper(self.hint, domains)
        D = daeVariable_wrapper(self.D, domains)
        Resistance = daeVariable_wrapper(self.Resistance, domains)

        kwall = self.kwall()
        Do = self.Do()
        Di = self.Di()
        pi = self.pi

        Resext = 1 / (pi * Do * hext)
        Resint = 1 / (pi * D * hint)
        Reswall = Log(Do / Di) / (2 * pi * kwall)

        eq.Residual = Resistance - (Resint + Reswall + Resext)


    def eq_total_he(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        Tsat = self.Tsat()
        T = daeVariable_wrapper(self.T, domains)
        Resistance = daeVariable_wrapper(self.Resistance, domains)
        Qout = daeVariable_wrapper(self.Qout, domains)

        eq.Residual = Qout * Resistance - (T - Tsat)


    def eq_calculate_To(self):

        eq = self.CreateEquation("WallHeat1", "Heat balance - wall")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        Do = self.Do()
        Tsat = self.Tsat()
        pi = self.pi

        To = daeVariable_wrapper(self.To, domains)
        hext = daeVariable_wrapper(self.hext, domains)
        Qout = daeVariable_wrapper(self.Qout, domains)

        eq.Residual = Qout - (To - Tsat) * (pi * Do * hext)


    def eq_calculate_Ti(self):

        eq = self.CreateEquation("WallHeat2", "Heat balance - wall")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        kwall = self.kwall()
        Do = self.Do()
        Di = self.Di()
        pi = self.pi

        Ti = daeVariable_wrapper(self.Ti, domains)
        To = daeVariable_wrapper(self.To, domains)
        Qout = daeVariable_wrapper(self.Qout, domains)

        eq.Residual = Qout *  Log(Do / Di) - (Ti - To) * (2 * pi * kwall)


    def DeclareEquations(self):

        daeModelExtended.DeclareEquations(self)

        self.eq_calculate_To()
        self.eq_calculate_Ti()
        self.eq_calculate_hint()
        self.eq_calculate_hext()
        self.eq_calculate_resistance()
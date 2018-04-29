import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, mol, J, W, rad


class Pipe(daeModel):

    def __init__(self, Name, Parent=None, Description=""):


        daeModel.__init__(self, Name, Parent, Description)

        # Defining Fluid
        self.fluid = daeThermoPackage("TPP", self, "")
        self.fluid.LoadCoolProp(["Water",],  # compound IDs in the mixture
            [],                              # compund CAS numbers (optional)
            {'Liquid': eLiquid},             # dictionary {'phaseLabel' : stateOfAggregation}
            eMass,                           # default basis is eMole (other options are eMass or eUndefinedBasis)
            {})                              # options dictionary (defaut is empty)
        self.fluid.fraction = [1.0,]

        # Domains

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Constants

        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa)

        # Defining Parameters

        self.tetha = daeParameter("tetha", rad, self, "Angle", [self.x, ])
        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")

        # Defining Variable Types

        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 373.0, 300.0, 0.01)
        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-10, 1e+02, 0.1, 1e-05)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        specific_heat_capacity_t = daeVariableType("specific_heat_capacity_t", (J ** (1))*(K ** (-1))*(kg ** (-1)), 1e-10, 1e+10, 1000., 1e-05)

        mass_flowrate_per_length_t = daeVariableType("mass_flowrate_per_length_t", (kg ** (1)) * (m ** (-1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)


        # Defining Variables

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", [self.x, ])
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.P = daeVariable("P", dimentionless_pressure_t, self, "Fluid Pressure", [self.x, ])
        self.v = daeVariable("v", velocity_t, self, "Fluid Velocity", [self.x, ])
        self.fD = daeVariable("f_D", darcy_t, self, "Darcy friction factor", [self.x, ])
        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Fluid thermal conductivity", [self.x, ])
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid", [self.x, ])
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid", [self.x, ])
        self.cp = daeVariable("h", specific_heat_capacity_t, self, "Heat capacity of the liquid", [self.x, ])

        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", [self.x, ])
        self.Mout = daeVariable("Mout", mass_flowrate_per_length_t, self, "Heat loss per length", [self.x, ])


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        # Definitions that simplify the notation
        dk_dx = lambda x: d(self.k(x), self.x, eCFDM)
        dT_dx = lambda x: d(self.T(x), self.x, eCFDM)

        # Equations: Fluid Property

        eq = self.CreateEquation("cp", "Heat capacity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.cp(x) - self.fluid.cp(self.PR*self.P(x),self.T(x),self.fluid.fraction)

        eq = self.CreateEquation("rho", "Density calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.rho(x) - self.fluid.rho(self.PR*self.P(x),self.T(x),self.fluid.fraction)

        eq = self.CreateEquation("mu", "Viscosity calculation.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.mu(x) - self.fluid.mu(self.PR*self.P(x),self.T(x),self.fluid.fraction)

        eq = self.CreateEquation("kappa", "Thermal conductivity calculation.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.kappa(x) - self.fluid.kappa(self.PR*self.P(x),self.T(x),self.fluid.fraction)

        # Equations: Fluidynamics

        eq = self.CreateEquation("v", "Velocity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D() ** 2
        eq.Residual = self.v(x) - self.k(x) / ( self.rho(x) * A )

        eq = self.CreateEquation("f_D", "Friction Factor")
        # TODO - Include laminar flow
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Re = self.D() * self.v(x) * self.rho(x) / self.mu(x)
        eq.Residual = 1 / Sqrt(self.fD(x)) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Sqrt(self.fD(x)))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * self.v(x) ** 2 / ( self.D() * self.g )
        DeltaP = self.g * self.rho(x) * hL
        # eq.Residual = dt(A*self.rho(x)*self.v(x)) + \
        #               d(self.k(x)*self.v(x), self.x, eCFDM)/self.L()  + \
        #               self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
        #               A*DeltaP + \
        #               A*self.rho(x)*self.g*Sin(self.tetha(x))
        eq.Residual = dt(A*self.rho(x)*self.v(x)) + \
                      self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
                      self.k(x) * d(self.v(x), self.x, eCFDM) / self.L() + \
                      A*DeltaP + \
                        A*self.rho(x)*self.g*Sin(self.tetha(x))

        eq = self.CreateEquation("MassBal", "Mass balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        eq.Residual = dt(A*self.rho(x)) + dk_dx(x) / self.L() + self.Mout(x)

        eq = self.CreateEquation("HeatBal", "Heat balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + d(self.k(x) * self.cp(x) * self.T(x), self.x,
        #                                                               eCFDM) / self.L()
        eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)


        # TODO Other equations

        # Para calcular Nusselt
        r"""Calculates the Nussel dimensionless number using Petukhov correlation modified by Gnielinski.
            See Incropera 4th Edition [8.63]
            :param \**kwargs:
                See below
            :Keyword Arguments:
                * *Re* (``float``) --
                    Reynolds number
                * *Pr* (``float``) --
                    Prandtl number
                * *f* (``float``) --
                    Fanning friction factor
            """
        #NuD = (f / 2.) * (Re - 1000.) * Pr / (1. + 12.7 * np.sqrt(f / 2.) * (np.power(Pr, 2 / 3) - 1.))

        # Para Calcular o Prandtl
        r"""Calculates the Prandtl dimensionless number
           :param \**kwargs:
               See below
           :Keyword Arguments:
               * *cp* (``float``) --
                   Heat capacity in kJ/kg
               * *mu* (``float``) --
                   Dynamic viscosity in Pa.s
               * *k* (``float``) --
                   Thermal conductivity in W/(m.K)
           """
        # Pr = 1000 * cp * mu / k

        # Para Calcular h_coefficient_internal_flow
        r"""Calculates heat transfer coefficient for internal flow
            in circular pipes in W/(m²K).
           :param \**kwargs:
               See below
           :Keyword Arguments:
               * *k* (``float``) --
                   Thermal conductivity in W/(m.K)
               * *D* (``float``) --
                   Pipe external diameter
               * *Nu* (``float``) --
                   Reynolds number
               """
        # h = NuD * k / D

        # Para calcular
        # def h_coefficient_external_condensation_unique_tube(**kwargs):
        #     r"""Calculates heat transfer coefficient for external condensation
        #      around an unique circular pipes in W/(m²K). See Incropera 4th Edition [10.41]
        #     :param \**kwargs:
        #         See below
        #     :Keyword Arguments:
        #         * *Ts* (``float``) --
        #             External pipe surface temperature in K
        #         * *Tsat* (``float``) --
        #             External saturated vapor temperature in K
        #         * *rho_l* (``float``) --
        #             Saturated liquid density in kg/m³
        #         * *mu_l* (``float``) --
        #             Liquid dynamic viscosity in Pa.s
        #         * *kl* (``float``) --
        #             Liquid thermal conductivity in W/(m.K)
        #         * *rho_v* (``float``) --
        #             Saturated vapor density in kg/m³
        #         * *hvap* (``float``) --
        #             Vaporization heat in J/kg
        #         * *D* (``float``) --
        #             Pipe external diameter
        #     """
        #
        #     g = const.physical_constants['standard acceleration of gravity'][0]
        #
        #     Ts = kwargs.pop('Ts')
        #     Tsat = kwargs.pop('Tsat')
        #     rho_l = kwargs.pop('rho_l')
        #     mu_l = kwargs.pop('mu_l')
        #     kl = kwargs.pop('kl')
        #     rho_v = kwargs.pop('rho_v')
        #     hvap = kwargs.pop('hvap')
        #     D = kwargs.pop('D')
        #
        #     num = (g * rho_l * (rho_l - rho_v) * kl ** 3. * hvap)
        #     den = mu_l * abs(Tsat - Ts) * D
        #     hd1 = 0.729 * (num / den) ** 0.25
        #
        #     # hd2 = 100 + Ts*0
        #     # hd = (hd1<hd2)*hd2 + (hd1>=hd2)*hd1
        #
        #     return hd1
        #
        # def h_coefficient_external_condensation(**kwargs):
        #     r"""Calculates heat transfer coefficient for external condensation
        #      around circular aligned pipes in W/(m²K). See Incropera 4th Edition [10.41]
        #     :param \**kwargs:
        #         See below
        #     :Keyword Arguments:
        #         * *Ts* (``float``) --
        #             External pipe surface temperature in K
        #         * *Tsat* (``float``) --
        #             External saturated vapor temperature in K
        #         * *rho_l* (``float``) --
        #             Saturated liquid density in kg/m³
        #         * *mu_l* (``float``) --
        #             Liquid dynamic viscosity in Pa.s
        #         * *kl* (``float``) --
        #             Liquid thermal conductivity in W/(m.K)
        #         * *rho_v* (``float``) --
        #             Saturated vapor density in kg/m³
        #         * *hvap* (``float``) --
        #             Vaporization heat in J/kg
        #         * *D* (``float``) --
        #             Pipe external diameter
        #         * *N* (``float``) --
        #             Number of pipes over the actual pipe, including it
        #     """
        #
        #     N = kwargs.pop('N')
        #
        #     hd1 = h_coefficient_external_condensation_unique_tube(**kwargs)
        #
        #     f = N ** 0.75 - (N - 1) ** 0.75
        #
        # return f * hd1

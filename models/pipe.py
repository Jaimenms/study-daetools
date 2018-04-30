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

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa)

        # Defining Parameters

        self.Text = daeParameter("T_ext", K, self, "External temperature")
        self.tetha = daeParameter("tetha", rad, self, "Angle", [self.x, ])
        self.D = daeParameter("D", m, self, "External pipe diameter")
        self.Dout = daeParameter("D_out", m, self, "External pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.kappawall = daeParameter("kappa_wall", (K ** (-1))*(W ** (1))*(m ** (-1)), self, "Wall thermal conductivity")

        # Defining Variable Types

        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-3, 1e+02, 0.1, 1e-05)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        specific_heat_capacity_t = daeVariableType("specific_heat_capacity_t", (J ** (1))*(K ** (-1))*(kg ** (-1)), 1e-10, 1e+10, 1000., 1e-05)

        mass_flowrate_per_length_t = daeVariableType("mass_flowrate_per_length_t", (kg ** (1)) * (m ** (-1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e-10, 1e+10, 0.1, 1e-05)

        # Defining Variables

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", [self.x, ])

        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.Tw1 = daeVariable("Tw_1", water_temperature_t, self, "Internal wall temperature", [self.x, ])
        self.Tw2 = daeVariable("Tw_2", water_temperature_t, self, "External wall temperature", [self.x, ])

        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Fluid thermal conductivity", [self.x, ])
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid", [self.x, ])
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid", [self.x, ])
        self.cp = daeVariable("h", specific_heat_capacity_t, self, "Heat capacity of the liquid", [self.x, ])

        self.v = daeVariable("v", velocity_t, self, "Fluid Velocity", [self.x, ])
        self.Re = daeVariable("Re", no_t, self, "Reynolds", [self.x, ])
        self.fD = daeVariable("f_D", darcy_t, self, "Darcy friction factor", [self.x, ])

        self.P = daeVariable("P", dimentionless_pressure_t, self, "Fluid Pressure", [self.x, ])

        self.prandtl = daeVariable("prandtl", no_t, self, "Prandtl", [self.x, ])
        self.nusselt = daeVariable("nusselt", no_t, self, "Nusselt", [self.x, ])
        self.hint = daeVariable("hint", heat_transfer_coefficient_t, self, "hint", [self.x, ])
        self.hext = daeVariable("hext", heat_transfer_coefficient_t, self, "hext", [self.x, ])

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
        eq.Residual = 1 / Pow(self.fD(x), 0.5) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Pow(self.fD(x),0.5))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * Pow(self.v(x), 2) / ( self.D() * self.g )
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

        # TODO Other equations

        eq = self.CreateEquation("Re", "Re")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Re(x) - self.D() * self.v(x) * self.rho(x) / self.mu(x)

        eq = self.CreateEquation("Prandtl", "Prandtl")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.prandtl(x) - self.cp(x) * self.mu(x) / self.kappa(x)

        eq = self.CreateEquation("Nusselt", "Nusselt")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.nusselt(x) - 0.027 * Pow(self.Re(x), 4/5) * Pow(self.prandtl(x), 1/3)
        # eq.Residual = self.nusselt(x) - (self.fD(x) / 8.) * (Re - 1000.) * self.prandtl(x) / (1. + 12.7 * Sqrt(self.fD(x) / 8.) * Pow(self.prandtl(x), 2 / 3) - 1.)

        eq = self.CreateEquation("hint", "hint")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.hint(x) - self.nusselt(x) * self.kappa(x) / self.D()

        eq = self.CreateEquation("HeatBal0", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + d(self.k(x) * self.cp(x) * self.T(x), self.x,
        #                                                              eCFDM) / self.L()
        eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L()
        #eq.Residual = d( self.T(x), self.x, eCFDM)


        eq = self.CreateEquation("HeatBal1", "Heat balance - Tw1")
        r"""Calculates the Nussel dimensionless number using Petukhov correlation modified by Gnielinski.
            See Incropera 4th Edition [8.63]
            """
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) - 2 * self.pi * self.D() * self.hint(x) * (self.T(x) - self.Tw1(x))

        eq = self.CreateEquation("HeatBal2", "Heat balance - Tw2")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) - 2 * self.pi * self.kappawall() * (self.Tw1(x) - self.Tw2(x))  / Log(self.Dout()/self.D())

        eq = self.CreateEquation("HeatBal3", "Heat balance - Q")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Qout(x) - 2 * self.pi * self.Dout() * self.hext(x) * (self.Tw2(x) - self.Text())

        eq = self.CreateEquation("ThermalRes", "Overall Thermal Resistance")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Rtotal  = 1/(2 * self.pi * self.Dout() * self.hext(x)) + \
                  1/(2 * self.pi * self.D() * self.hint(x)) + \
                  Log(self.Dout()/self.D()) / (2 * self.pi * self.kappawall())
        eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / Rtotal


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

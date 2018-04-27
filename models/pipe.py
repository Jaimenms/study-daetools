import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, mol, J, W


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

        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")

        # Defining Variable Types

        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 373.0, 300.0, 0.01)
        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-10, 1e+02, 0.1, 1e-05)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)

        # Defining Variables

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", [self.x, ])
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.P = daeVariable("P", dimentionless_pressure_t, self, "Fluid Pressure", [self.x, ])
        self.v = daeVariable("v", velocity_t, self, "Fluid Velocity", [self.x, ])
        self.fD = daeVariable("f_D", darcy_t, self, "Darcy friction factor", [self.x, ])
        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Fluid thermal conductivity", [self.x, ])
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid", [self.x, ])
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid", [self.x, ])

    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        # Definitions that simplify the notation
        dP_dx = lambda x: d(self.P(x), self.x, eCFDM)
        dk_dx = lambda x: d(self.k(x), self.x, eCFDM)
        dT_dx = lambda x: d(self.T(x), self.x, eCFDM)

        # Equations: Fluid Property

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
        DeltaP = self.g * self.rho(x) * self.L() * hL / self.PR
        eq.Residual = dP_dx(x) + DeltaP

        eq = self.CreateEquation("MassBal", "Mass balance")
        # TODO - Include balances
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        eq.Residual = dk_dx(x)

        eq = self.CreateEquation("HeatBal", "Heat balance")
        # TODO - Include balances
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        eq.Residual = dT_dx(x)

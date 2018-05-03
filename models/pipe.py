import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, mol, J, W, rad


class ConvectionPipe(daeModel):

    def __init__(self, Name, Parent=None, Description=""):


        daeModel.__init__(self, Name, Parent, Description)

        # Domains

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Constants

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa)

        # Defining Parameters

        self.tetha = daeParameter("tetha", rad, self, "Angle", [self.x, ])
        self.D = daeParameter("D", m, self, "Internal pipe diameter")
        self.Dout = daeParameter("Dout", m, self, "External pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.kwall = daeParameter("k_wall", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-1)), self, "Wall conductivity")
        self.Text = daeParameter("T_ext", K, self, "External Temperature")
        self.hext = daeParameter("h_ext", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-2)), self, "External heat transfer coefficient")


        # Defining Variable Types

        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10,
                                            0.1, 1e-05)

        # Defining Variables

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", [self.x, ])
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.Tw = daeVariable("Tw", water_temperature_t, self, "Wall Temperature", [self.x, ])
        self.P = daeVariable("P", dimentionless_pressure_t, self, "Fluid Pressure", [self.x, ])

        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid")
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid")
        self.cp = daeVariable("cp", specific_heat_capacity_t, self, "Heat capacity of the liquid")
        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Thermal Conductivity of the liquid")

        self.fD = daeVariable("f_D", darcy_t, self, "Darcy friction factor", [self.x, ])
        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", [self.x, ])


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        # Definitions that simplify the notation
        dk_dx = lambda x: d(self.k(x), self.x, eCFDM)

        # Equations: Fluidynamics


        eq = self.CreateEquation("f_D", "Friction Factor")
        # TODO - Include laminar flow
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D() ** 2
        v = self.k(x) / ( self.rho() * A )
        Re = self.D() * v * self.rho() / self.mu()
        eq.Residual = 1 / Pow(self.fD(x), 0.5) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Pow(self.fD(x),0.5))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * Pow(v, 2) / ( self.D() * self.g )
        DeltaP = self.g * self.rho() * hL
        # eq.Residual = dt(A*self.rho(x)*self.v(x)) + \
        #               d(self.k(x)*self.v(x), self.x, eCFDM)/self.L()  + \
        #               self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
        #               A*DeltaP + \
        #               A*self.rho(x)*self.g*Sin(self.tetha(x))
        eq.Residual = self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
                      A*DeltaP + \
                      A*self.rho()*self.g*Sin(self.tetha(x))

        eq = self.CreateEquation("MassBal", "Mass balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        # eq.Residual = dt(A*self.rho(x)) + dk_dx(x) / self.L() + self.Mout(x)
        eq.Residual = dk_dx(x)

        eq = self.CreateEquation("HeatBal0", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)
        eq.Residual = A * self.rho() * self.cp() * dt(self.T(x)) + self.k(x) * self.cp() * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)

        eq = self.CreateEquation("IntConv", "IntConv")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        prandtl = self.cp() * self.mu() / self.kappa()
        nusselt =  0.027 * (Re ** (4/5)) * (prandtl ** (1/3))
        #nusselt = (self.fD(x) / 8.) * (Re - 1000.) * prandtl / (
        #            1. + 12.7 * Sqrt(self.fD(x) / 8.) * Pow(prandtl, 2 / 3) - 1.)
        hint = nusselt * self.kappa() / self.D()
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.D() * hext)
        Resint = 1 / (2 * self.pi * self.Dout() * hint)
        Reswall = (self.D() - self.Dout()) / ( 4 * self.pi * self.Dout() * self.kwall() )
        eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / (Resint + Reswall + Resext)
        #eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / (Reswall)

        eq = self.CreateEquation("ExtConv", "ExtConv")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.Dout() * hext)
        eq.Residual = self.Qout(x) - (self.Tw(x) - self.Text()) / Resext



class IsothermalPipe(daeModel):

    def __init__(self, Name, Parent=None, Description=""):


        daeModel.__init__(self, Name, Parent, Description)

        # Domains

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Constants

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa)

        # Defining Parameters

        self.tetha = daeParameter("tetha", rad, self, "Angle", [self.x, ])
        self.D = daeParameter("D", m, self, "External pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")

        # Defining Variable Types

        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)

        # Defining Variables

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", [self.x, ])
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.P = daeVariable("P", dimentionless_pressure_t, self, "Fluid Pressure", [self.x, ])

        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid")
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid")

        self.fD = daeVariable("f_D", darcy_t, self, "Darcy friction factor", [self.x, ])


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        # Definitions that simplify the notation
        dk_dx = lambda x: d(self.k(x), self.x, eCFDM)
        dT_dx = lambda x: d(self.T(x), self.x, eCFDM)

        # Equations: Fluidynamics


        eq = self.CreateEquation("f_D", "Friction Factor")
        # TODO - Include laminar flow
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D() ** 2
        v = self.k(x) / ( self.rho() * A )
        Re = self.D() * v * self.rho() / self.mu()
        eq.Residual = 1 / Pow(self.fD(x), 0.5) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Pow(self.fD(x),0.5))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * Pow(v, 2) / ( self.D() * self.g )
        DeltaP = self.g * self.rho() * hL
        # eq.Residual = dt(A*self.rho(x)*self.v(x)) + \
        #               d(self.k(x)*self.v(x), self.x, eCFDM)/self.L()  + \
        #               self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
        #               A*DeltaP + \
        #               A*self.rho(x)*self.g*Sin(self.tetha(x))
        eq.Residual = self.PR*d(A*self.P(x), self.x, eCFDM)/self.L() + \
                      A*DeltaP + \
                      A*self.rho()*self.g*Sin(self.tetha(x))

        eq = self.CreateEquation("MassBal", "Mass balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        # eq.Residual = dt(A*self.rho(x)) + dk_dx(x) / self.L() + self.Mout(x)
        eq.Residual = dk_dx(x)

        eq = self.CreateEquation("HeatBal0", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)
        eq.Residual = dT_dx(x)
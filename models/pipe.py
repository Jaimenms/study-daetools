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

        # Domains

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Constants

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa)

        # Defining Parameters

        self.tetha = daeParameter("tetha", rad, self, "Angle", [self.x, ])
        self.Di = daeParameter("D_i", m, self, "Inside pipe diameter")
        self.Do = daeParameter("D_o", m, self, "Outside pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.kwall = daeParameter("k_wall", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-1)), self, "Wall conductivity")
        self.Text = daeParameter("T_ext", K, self, "External Temperature")
        self.hext = daeParameter("h_ext", (K ** (-1))*(J ** (1))*(s ** (-1))*(m ** (-2)), self, "External heat transfer coefficient")


        # Defining Variable Types

        diameter_t = daeVariableType("diameter_t", (m ** (1)), 0.001, 1.0, 0.5, 1e-05)
        dimentionless_pressure_t = daeVariableType("dimentionless_pressure_t", dimless, 0.5, 20.0, 1.0, 1e-06)
        darcy_t = daeVariableType("darcy_t", dimless, 1e-10, 100.0, 0.05, 1e-06)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-10, 1e+04, 0.1, 1e-05)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10,
                                            0.1, 1e-05)
        mass_biofilm_t = daeVariableType("mass_biofilm_t", (kg ** (1)) * (m ** (-2)), 1e-10, 1e+04, 0.1, 1e-05)

        # Defining Variables

        self.mf = daeVariable("m_f", mass_biofilm_t, self, "Mass of Biofilm", [self.x, ])
        self.D = daeVariable("D", diameter_t, self, "Internal flow diameter", [self.x, ])
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

        eq = self.CreateEquation("D", "Internal flow diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.D(x) - self.Di()


        eq = self.CreateEquation("f_D", "Friction Factor")
        # TODO - Include laminar flow
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D(x) ** 2
        v = self.k(x) / ( self.rho() * A )
        Re = self.D(x) * v * self.rho() / self.mu()
        eq.Residual = 1 / Pow(self.fD(x), 0.5) + 2. * Log10(self.ep() / 3.7 / self.D(x) + 2.51 / Re / Pow(self.fD(x),0.5))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * Pow(v, 2) / ( self.D(x) * self.g )
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

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = dt(A * self.rho(x) * self.cp(x) * self.T(x)) + self.k(x) * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)
        eq.Residual = A * self.rho() * self.cp() * dt(self.T(x)) + self.k(x) * self.cp() * d( self.T(x), self.x, eCFDM) / self.L() + self.Qout(x)

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        prandtl = self.cp() * self.mu() / self.kappa()
        nusselt =  0.027 * (Re ** (4/5)) * (prandtl ** (1/3))
        #nusselt = (self.fD(x) / 8.) * (Re - 1000.) * prandtl / (
        #            1. + 12.7 * Sqrt(self.fD(x) / 8.) * Pow(prandtl, 2 / 3) - 1.)
        hint = nusselt * self.kappa() / self.D(x)
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.Di() * hext)
        Resint = 1 / (2 * self.pi * self.D(x) * hint)
        Reswall = Log(self.Do()/self.Di()) / (2 * self.pi * self.kwall() )
        eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / (Resint + Reswall + Resext)

        eq = self.CreateEquation("WallHeat", "Heat balance - wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        hext = self.hext()
        Resext = 1 / (2 * self.pi * self.Do() * hext)
        eq.Residual = self.Qout(x) - (self.Tw(x) - self.Text()) / Resext


        self.IF(Time() < Constant(5*s), eventTolerance = 1E-5)

        eq = self.CreateEquation("Biofilm", "Biofilm Formation - OFF")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.mf(x)

        self.ELSE()

        eq = self.CreateEquation("Biofilm", "Biofilm Formation - ON")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        vast = v / Constant(1 * m / s)
        Tast = self.T(x) / Constant(1 * K)
        Jp = (1.19e-7 - 1.14e-7 * vast) * Constant( 1 * kg / m **2 / s)
        b =  1 / (4.26e4  + 3.16e5 * vast) * Constant( 1 / s)
        Jr = b * self.mf(x)
        # k27_mf = 0.599
        # sigmoid = 1 * Exp(0.6221 * (Tast - 315.34)) / (1 + Exp(0.6221 * (Tast - 315.34)))
        # k_mf = 1.624e7 * Tast * Exp(-1000 * 13.609 / 1.987 / Tast) * (1 - sigmoid)
        eq.Residual = dt(self.mf(x)) - ( Jp - Jr )
        # eq.Residual = dt(self.mf(x)) - ( Jp * k27_mf / k_mf * Constant(1. * kg * m ** (-2) * s ** (-1) ) - Jr * Constant(1. * s ** (-1) ))

        self.END_IF()


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
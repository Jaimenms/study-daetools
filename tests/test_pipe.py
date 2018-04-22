import os
import sys
import inspect
import json
import pprint

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from daetools_reporters import setupDataReporters


class Pipe(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        from pyUnits import m, kg, s, K, Pa, mol, J, W

        daeModel.__init__(self, Name, Parent, Description)

        self.x = daeDomain("x", self, unit(), "Axial Position")

        # Defining Parameters

        self.g = daeParameter("g", m / (s **2), self, "Gravity acceleration")
        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.rho = daeParameter("rho", kg / (m ** 3), self, "Density")
        self.mu = daeParameter("mu", Pa * s, self, "Dynamic viscosity")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.PR = daeParameter("PR", Pa, self, "Pressure reference")

        # self.PA = daeParameter("PA", unit(), self, "Pressure Inlet")
        # self.PB = daeParameter("PB", unit(), self, "Pressure Outlet")

        # Defining Variables
        self.P = daeVariable("P", no_t, self)
        self.P.Description = "Pressure"
        self.P.DistributeOnDomain(self.x)

        self.v = daeVariable("v", velocity_t, self)
        self.v.Description = "Velocity"

        self.fD = daeVariable("f_D", no_t, self)
        self.fD.Description = "Darcy friction factor"

        self.k = daeVariable("k", mass_flowrate_t, self)
        self.k.Description = "Mass flowrate"


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        # Intermediate variables
        A = 0.25 * 3.14 * self.D() ** 2
        Re = self.D() * self.v() * self.rho() / self.mu()
        hL = 0.5 * self.fD() * self.v() ** 2 / ( self.D() * self.g() )
        DeltaP = self.g() * self.rho() * hL * self.L() / self.PR()

        # Equations
        eq = self.CreateEquation("v", "Velocity")
        eq.Residual = self.v() - self.k() / ( self.rho() * A )

        eq = self.CreateEquation("f_D", "Friction Factor")
        eq.Residual = 1 / Sqrt(self.fD()) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Sqrt(self.fD()))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        eq.Residual = self.P(0) - self.P(1) - DeltaP


class sim_test1(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = Pipe("test1")
        self.m.Description = "Testing the solution for a pipe"
        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        self.m.x.CreateArray(2)

        # Setting Parameter values
        self.m.g.SetValue( 9.81 * m/ s**2 )
        self.m.D.SetValue( 4.026*0.0254 * m )
        self.m.L.SetValue( 1.0 * m )
        self.m.rho.SetValue( 1000 * kg / m**3 )
        self.m.mu.SetValue( 0.001 * Pa * s )
        self.m.ep.SetValue( 0.0018*0.0254 * m )
        self.m.PR.SetValue( 100000. * Pa )


    def SetUpVariables(self):

        # Setting Variable Initial Guesses
        self.m.fD.SetInitialGuesses(0.018 * unit())
        self.m.v.SetInitialGuesses(2. * m / s)

        self.m.k.SetInitialGuesses(10. * kg / s)
        self.m.P.AssignValue(0, 1. * unit())
        self.m.P.AssignValue(1, 0.995 * unit())


    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)


def test_test1():

    data = main()

    assert round( data['v']['Values'][0] - 2.36909, 4) == 0.0


def main():

    log = daeStdOutLog()
    solver = daeIDAS()

    #cfg = daeGetConfig()
    #cfg.SetBoolean("daetools.IDAS.printInfo", True)
    #cfg.SetBoolean("daetools.IDAS.LineSearchOffIC", True)
    #cfg.SetFloat("daetools.IDAS.MaxStep", 0.0)
    #cfg.SetFloat("daetools.IDAS.NonlinConvCoef", 0.1)
    #cfg.SetInteger("daetools.IDAS.numberOfSTNRebuildsDuringInitialization", 10000)
    #print(cfg)

    simulation = sim_test1()
    simulation.m.SetReportingOn(True)
    #simulation.m.SetReportingInterval(10)
    #simulation.m.SetTimeHorizon(1000)

    datareporter = setupDataReporters(simulation)

    simulation.Initialize(solver, datareporter, log)

    # Solve at time = 0 (initialization)
    simulation.SolveInitial()

    # Run
    simulation.Run()

    # Clean up
    simulation.Finalize()

    return json.load(open(simulation.report_filename))


if __name__ == "__main__":

    pprint.pprint(main(), indent = 4)
import os
import sys
import inspect
import json

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from daetools_reporters import setupDataReporters


class mod_test1(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        from pyUnits import m, kg, s, K, Pa, mol, J, W

        daeModel.__init__(self, Name, Parent, Description)

        self.g = daeParameter("g", m / (s **2), self, "Gravity acceleration")
        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.rho = daeParameter("rho", kg / (m ** 3), self, "Density")
        self.mu = daeParameter("mu", Pa * s, self, "Dynamic viscosity")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.k0 = daeParameter("k0", kg/s, self, "Mass flow ate")

        self.v = daeVariable("v", velocity_t, self)
        self.v.Description = "Velocity"

        self.Re = daeVariable("Re", no_t, self)
        self.Re.Description = "Reynolds"

        self.hL = daeVariable("h_L", no_t, self)
        self.hL.Description = "Head Loss"

        self.fD = daeVariable("f_D", no_t, self)
        self.fD.Description = "Darcy friction factor"

        self.deltaP = daeVariable("deltaP", pressure_t, self)
        self.deltaP.Description = "Pressure loss"

        self.k = daeVariable("k", mass_flowrate_t, self)
        self.k.Description = "Mass flowrate"


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        A = 0.25 * 3.14 * self.D() ** 2

        eq = self.CreateEquation("v", "Velocity")
        eq.Residual = self.v() - self.k() / ( self.rho() * A )

        eq = self.CreateEquation("Re", "Reynolds")
        eq.Residual = self.Re() - self.D() * self.v() * self.rho() / self.mu()

        eq = self.CreateEquation("f_D", "Friction Factor")
        eq.Residual = 1 / Sqrt(self.fD()) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / self.Re() / Sqrt(self.fD()))

        eq = self.CreateEquation("hL", "Head loss")
        eq.Residual = self.hL() - 0.5 * self.fD() * self.v() ** 2 / ( self.D() * self.g() )

        eq = self.CreateEquation("MomBal", "Momentum balance")
        eq.Residual = self.deltaP() - self.g() * self.rho() * self.hL() * self.L()

        eq = self.CreateEquation("MasBal", "Mass")
        eq.Residual = self.k() - self.k0()


class sim_test1(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = mod_test1("test1")
        self.m.Description = "Testing the solution for a pipe"
        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        self.m.g.SetValue(9.81*m/(s**2))
        self.m.D.SetValue(4.026*0.0254*m)
        self.m.L.SetValue(1.0*m)
        self.m.rho.SetValue(1000*kg/(m**3))
        self.m.mu.SetValue(0.001*Pa*s)
        self.m.ep.SetValue(0.0018*0.0254*m)
        self.m.k0.SetValue(22.*kg/s)


    def SetUpVariables(self):

        self.m.Re.SetInitialGuesses(285520. * unit())
        self.m.fD.SetInitialGuesses(0.001 * unit())
        self.m.hL.SetInitialGuesses(0.001 * unit())
        self.m.v.SetInitialGuesses(2. * m / s)
        self.m.deltaP.SetInitialGuesses(400. * Pa)

        #self.InitialConditionMode = eQuasiSteadyState


    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)


def test_test1():

    sim = sim_test1()

    daeActivity.simulate(
        sim,
        reportingInterval=5,
        timeHorizon=500,
        datareporter=setupDataReporters(sim),
        guiRun=False)

    data = json.load(open(sim.report_filename))

    assert data['T']['Values'][0] == 310.0



if __name__ == "__main__":

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

    data = json.load(open(simulation.report_filename))

    print(data)
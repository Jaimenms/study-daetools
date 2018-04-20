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


        self.Q = daeParameter("En", J, self, "Heat")
        self.m = daeParameter("m", kg, self, "Mass")
        self.cp = daeParameter("c_p", J / (kg * K), self, "Specific heat capacity")
        self.T0 = daeParameter("T0", K, self, "Initial temperature")

        self.T = daeVariable("T", temperature_t, self)
        self.T.Description = "Final Temperature"


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("HeatBalance", "Heat balance")
        eq.Residual = self.Q() - self.m() * self.cp() * ( self.T() - self.T0())


class sim_test1(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = mod_test1("test1")
        self.m.Description = "Testing the solution of a linear equation"
        self.report_filename = __file__ + '.json'

    def SetUpParametersAndDomains(self):

        self.m.Q.SetValue(10 * 385 * J)
        self.m.cp.SetValue(385 * J / (kg * K))
        self.m.m.SetValue(1 * kg)
        self.m.T0.SetValue(300 * K)


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

    simulation = sim_test1()

    daeActivity.simulate(simulation, reportingInterval=5, timeHorizon=500,
                         datareporter=setupDataReporters(simulation), guiRun=False)

    data = json.load(open(simulation.report_filename))

    print(data)
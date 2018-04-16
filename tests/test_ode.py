import os
import sys
import inspect
import json

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *

from daetools_reporters import setupDataReporters


class mod_test2(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        from pyUnits import m, kg, s, K, Pa, mol, J, W

        daeModel.__init__(self, Name, Parent, Description)


        self.Q = daeParameter("Q", W, self, "Power")
        self.m = daeParameter("m", kg, self, "Mass")
        self.cp = daeParameter("c_p", J / (kg * K), self, "Specific heat capacity")

        self.T = daeVariable("T", temperature_t, self)
        self.T.Description = "Temperature"


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("HeatBalance", "Heat balance")
        eq.Residual = dt(self.m() * self.cp() * self.T()) - self.Q()



class sim_test2(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = mod_test2("test1")
        self.m.Description = "Testing the solution of a linear equation"
        self.report_filename = __file__ + '.json'

    def SetUpParametersAndDomains(self):

        self.m.Q.SetValue(385 * W)
        self.m.cp.SetValue(385 * J / (kg * K))
        self.m.m.SetValue(1 * kg)

    def SetUpVariables(self):
        # SetInitialCondition function in the case of distributed variables can accept additional arguments
        # specifying the indexes in the domains. In this example we loop over the open x and y domains,
        # thus we start the loop with 1 and end with NumberOfPoints-1 (for both domains)
        self.m.T.SetInitialCondition(300 * K)


def test_test2():

    sim = sim_test2()

    daeActivity.simulate(
        sim,
        reportingInterval=1,
        timeHorizon=10,
        datareporter=setupDataReporters(sim),
        guiRun=False)

    data = json.load(open(sim.report_filename))

    from scipy.interpolate import interp1d
    from numpy import asscalar
    f = interp1d(data['T']['Times'],data['T']['Values'])

    assert round(asscalar(f(10)) - 310.0, 4) == 0.0



if __name__ == "__main__":

    simulation = sim_test2()

    daeActivity.simulate(simulation, reportingInterval=1, timeHorizon=10,
                         datareporter=setupDataReporters(simulation), guiRun=False)

    data = json.load(open(simulation.report_filename))

    print(data)
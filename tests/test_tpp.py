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

from pyUnits import m, kg, s, K, Pa, mol, J, W


class PipeWithTPP(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        daeModel.__init__(self, Name, Parent, Description)

        self.tpp = daeThermoPackage("TPP", self, "")

        self.tpp.LoadCoolProp(["Water",],       # compound IDs in the mixture
            [],                                 # compund CAS numbers (optional)
            {'Liquid': eLiquid},                # dictionary {'phaseLabel' : stateOfAggregation}
            eMass,                              # default basis is eMole (other options are eMass or eUndefinedBasis)
            {"backend": "HEOS",
            "referenceState": "DEF",
            "debugLevel": "0"})                                 # options dictionary (defaut is empty)


        self.T = daeParameter("T", K, self, "Temperature")
        self.rho   = daeVariable("rho",     density_t,                self, "Density of the liquid")




    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        #dP_dx = lambda x: d(self.P(x), self.x, eCFDM)

        P_aux = 1e5 * Pa
        x_aux = [1.0,]
        T_aux = 300 * K

        # Equations
        eq = self.CreateEquation("rho", "Density calculation")
        #x = eq.DistributeOnDomain(self.x, eClosedClosed)
        #eq.Residual = self.rho(x) - self.tpp.rho(self.PR()*self.P(x),self.T(),[1.0,])
        eq.Residual = self.rho() - self.tpp.rho(P_aux, self.T(),x_aux)
        #eq.Residual = self.rho() - Constant(1000. * kg / m ** 3)


class sim_test1b(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = PipeWithTPP("test1b")
        self.m.Description = "Testing the solution for a pipe with Thermal Package"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):


        #self.m.x.CreateStructuredGrid(10, 0.0, 1.0)

        # Setting Parameter values
        self.m.T.SetValue( 300 * K)


    def SetUpVariables(self):

        # Setting Variable Initial Guesses
        self.m.rho.SetInitialGuesses(1000. * kg / m ** 3)

        #Nx = self.m.x.NumberOfPoints
        #for i in range(Nx):
        #    self.m.P.AssignValue(i, 1. * unit())

        pass


    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)




def main(simulation = None):

    log = daeStdOutLog()
    solver = daeIDAS()

    cfg = daeGetConfig()
    cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')

    simulation.m.SetReportingOn(True)

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

    pprint.pprint(main(simulation = sim_test1b()), indent = 4)
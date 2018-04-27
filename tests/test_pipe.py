import os
import sys
import inspect
import json
import pprint
import argparse

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from daetools_reporters import setupDataReporters


class Pipe(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        from pyUnits import m, kg, s, K, Pa, mol, J, W

        daeModel.__init__(self, Name, Parent, Description)

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Parameters

        #self.g = daeParameter("g", m / (s **2), self, "Gravity acceleration")
        self.g = Constant(9.81 * m / s **2 )
        self.PR = Constant(1e5 * Pa)
        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.rho = daeParameter("rho", kg / (m ** 3), self, "Density")
        self.mu = daeParameter("mu", Pa * s, self, "Dynamic viscosity")
        self.ep = daeParameter("epsilon", m, self, "Roughness")

        # Defining Variables
        self.P = daeVariable("P", no_t, self, description = "Pressure", domains=[self.x,])

        self.v = daeVariable("v", velocity_t, self)
        self.v.Description = "Velocity"

        self.fD = daeVariable("f_D", no_t, self)
        self.fD.Description = "Darcy friction factor"

        self.k = daeVariable("k", mass_flowrate_t, self)
        self.k.Description = "Mass flowrate"


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        dP_dx = lambda x: d(self.P(x), self.x, eCFDM)

        # Intermediate variables
        A = 0.25 * 3.14 * self.D() ** 2
        Re = self.D() * self.v() * self.rho() / self.mu()
        hL = 0.5 * self.fD() * self.v() ** 2 / ( self.D() * self.g )
        DeltaP = self.g * self.rho() * self.L() * hL / self.PR

        # Equations
        eq = self.CreateEquation("v", "Velocity")
        eq.Residual = self.v() - self.k() / ( self.rho() * A )

        eq = self.CreateEquation("f_D", "Friction Factor")
        eq.Residual = 1 / Sqrt(self.fD()) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Sqrt(self.fD()))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        #eq.Residual = self.P(0) - self.P(1) - DeltaP
        eq.Residual = dP_dx(x) + DeltaP



class TwoPipes(daeModel):

    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self, Name, Parent, Description)

        self.pipe = {}
        self.pipe['pipe_0'] = Pipe("pipe_0", self)
        self.pipe['pipe_1'] = Pipe("pipe_1", self)

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)


class TwoConnectedPipes(daeModel):

    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self, Name, Parent, Description)

        # Defining Variables
        self.Pnode = daeVariable("P_node", no_t, self)
        self.Pnode.Description = "Nodal Pressure"

        self.pipe = []
        self.pipe.append(Pipe("pipe_0", self))
        self.pipe.append(Pipe("pipe_1", self))



    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("Connection_Node_1_a")
        eq.Residual = self.Pnode() - self.pipe[1].P(0)

        eq = self.CreateEquation("Connection_Node_1_b")
        eq.Residual = self.pipe[0].P(10) - self.Pnode()

        eq = self.CreateEquation("Connection_Node_2")
        eq.Residual = self.pipe[0].k() - self.pipe[1].k()



class sim_test1a(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = Pipe("test1a")
        self.m.Description = "Testing the solution for a pipe"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        #self.m.x.CreateArray(2)
        self.m.x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        #self.m.g.SetValue( 9.81 * m/ s**2 )
        self.m.D.SetValue( 4.026*0.0254 * m )
        self.m.L.SetValue( 1.0 * m )
        self.m.rho.SetValue( 1000 * kg / m**3 )
        self.m.mu.SetValue( 0.001 * Pa * s )
        self.m.ep.SetValue( 0.0018*0.0254 * m )

    def SetUpVariables(self):

        # Setting Variable Initial Guesses
        self.m.fD.SetInitialGuesses(0.018 * unit())
        self.m.v.SetInitialGuesses(2. * m / s)
        self.m.k.SetInitialGuesses(10. * kg / s)
        self.m.P.AssignValue(0, 1. * unit())
        self.m.P.AssignValue(10, 0.995 * unit())


    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)


class sim_test2(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = TwoPipes("test2")
        self.m.Description = "Testing the solution for two disconnected pipes"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        #self.m.x.CreateArray(2)
        self.m.pipe['pipe_0'].x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        self.m.pipe['pipe_0'].g.SetValue( 9.81 * m/ s**2 )
        self.m.pipe['pipe_0'].D.SetValue( 4.026*0.0254 * m )
        self.m.pipe['pipe_0'].L.SetValue( 1.0 * m )
        self.m.pipe['pipe_0'].rho.SetValue( 1000 * kg / m**3 )
        self.m.pipe['pipe_0'].mu.SetValue( 0.001 * Pa * s )
        self.m.pipe['pipe_0'].ep.SetValue( 0.0018*0.0254 * m )
        self.m.pipe['pipe_0'].PR.SetValue( 100000. * Pa )

        #self.m.x.CreateArray(2)
        self.m.pipe['pipe_1'].x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        self.m.pipe['pipe_1'].g.SetValue( 9.81 * m/ s**2 )
        self.m.pipe['pipe_1'].D.SetValue( 4.026*0.0254 * m )
        self.m.pipe['pipe_1'].L.SetValue( 1.0 * m )
        self.m.pipe['pipe_1'].rho.SetValue( 1000 * kg / m**3 )
        self.m.pipe['pipe_1'].mu.SetValue( 0.001 * Pa * s )
        self.m.pipe['pipe_1'].ep.SetValue( 0.0018*0.0254 * m )
        self.m.pipe['pipe_1'].PR.SetValue( 100000. * Pa )

    def SetUpVariables(self):

        # Setting Variable Initial Guesses
        self.m.pipe['pipe_0'].fD.SetInitialGuesses(0.018 * unit())
        self.m.pipe['pipe_0'].v.SetInitialGuesses(2. * m / s)
        self.m.pipe['pipe_0'].k.SetInitialGuesses(10. * kg / s)
        self.m.pipe['pipe_0'].P.AssignValue(0, 1. * unit())
        self.m.pipe['pipe_0'].P.AssignValue(10, 0.995 * unit())

        # Setting Variable Initial Guesses
        self.m.pipe['pipe_1'].fD.SetInitialGuesses(0.018 * unit())
        self.m.pipe['pipe_1'].v.SetInitialGuesses(2. * m / s)
        self.m.pipe['pipe_1'].k.SetInitialGuesses(10. * kg / s)
        self.m.pipe['pipe_1'].P.AssignValue(0, 1. * unit())
        self.m.pipe['pipe_1'].P.AssignValue(10, 0.995 * unit())

    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)


class sim_test3(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = TwoConnectedPipes("test3")
        self.m.Description = "Testing the solution for two connected pipes"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        #self.m.x.CreateArray(2)
        self.m.pipe[0].x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        self.m.pipe[0].g.SetValue( 9.81 * m/ s**2 )
        self.m.pipe[0].D.SetValue( 4.026*0.0254 * m )
        self.m.pipe[0].L.SetValue( 1.0 * m )
        self.m.pipe[0].rho.SetValue( 1000 * kg / m**3 )
        self.m.pipe[0].mu.SetValue( 0.001 * Pa * s )
        self.m.pipe[0].ep.SetValue( 0.0018*0.0254 * m )
        self.m.pipe[0].PR.SetValue( 100000. * Pa )

        #self.m.x.CreateArray(2)
        self.m.pipe[1].x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        self.m.pipe[1].g.SetValue( 9.81 * m/ s**2 )
        self.m.pipe[1].D.SetValue( 4.026*0.0254 * m )
        self.m.pipe[1].L.SetValue( 1.0 * m )
        self.m.pipe[1].rho.SetValue( 1000 * kg / m**3 )
        self.m.pipe[1].mu.SetValue( 0.001 * Pa * s )
        self.m.pipe[1].ep.SetValue( 0.0018*0.0254 * m )
        self.m.pipe[1].PR.SetValue( 100000. * Pa )

    def SetUpVariables(self):

        self.m.Pnode.SetInitialGuesses(0.998 * unit())

        # Setting Variable Initial Guesses
        self.m.pipe[0].fD.SetInitialGuesses(0.018 * unit())
        self.m.pipe[0].v.SetInitialGuesses(2. * m / s)
        self.m.pipe[0].k.SetInitialGuesses(10. * kg / s)
        self.m.pipe[0].P.AssignValue(0, 1. * unit())

        # Setting Variable Initial Guesses
        self.m.pipe[1].fD.SetInitialGuesses(0.018 * unit())
        self.m.pipe[1].v.SetInitialGuesses(2. * m / s)
        self.m.pipe[1].k.SetInitialGuesses(10. * kg / s)
        self.m.pipe[1].P.AssignValue(10, 0.995 * unit())


    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)





class PipeWithTPP(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        from pyUnits import m, kg, s, K, Pa, mol, J, W

        daeModel.__init__(self, Name, Parent, Description)

        self.tpp = daeThermoPackage("TPP", self, "")

        self.x = daeDomain("x", self, unit(), "X axis domain")

        # Defining Parameters

        self.T = daeParameter("T", K, self, "Temperature")
        self.g = daeParameter("g", m / (s **2), self, "Gravity acceleration")
        self.D = daeParameter("D", m, self, "Diameter")
        self.L = daeParameter("L", m, self, "Length")
        #self.rho = daeParameter("rho", kg / (m ** 3), self, "Density")
        #self.mu = daeParameter("mu", Pa * s, self, "Dynamic viscosity")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.PR = daeParameter("PR", Pa, self, "Pressure reference")

        # Defining Variables
        self.P = daeVariable("P", no_t, self, "Pressure")
        self.P.DistributeOnDomain(self.x)

        self.v = daeVariable("v", velocity_t, self, "Velocity")
        self.v.DistributeOnDomain(self.x)

        self.fD = daeVariable("f_D", no_t, self, "Darcy friction factor")
        self.fD.DistributeOnDomain(self.x)

        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate")

        self.kappa = daeVariable("kappa",     thermal_conductivity_t,   self, "Thermal conductivity of the liquid")
        self.kappa.DistributeOnDomain(self.x)

        self.mu    = daeVariable("mu",    dynamic_viscosity_t,      self, "Viscosity of the liquid")
        self.mu.DistributeOnDomain(self.x)

        self.rho   = daeVariable("rho",     density_t,                self, "Density of the liquid")
        self.rho.DistributeOnDomain(self.x)

        self.tpp.LoadCoolProp(["Water",],       # compound IDs in the mixture
            [],                                 # compund CAS numbers (optional)
            {'Liquid': eLiquid},                # dictionary {'phaseLabel' : stateOfAggregation}
            eMass,                              # default basis is eMole (other options are eMass or eUndefinedBasis)
            {"backend": "HEOS",
            "referenceState": "DEF",
            "debugLevel": "0"})                                 # options dictionary (defaut is empty)



    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        dP_dx = lambda x: d(self.P(x), self.x, eCFDM)

        ctpp = daeThermoPackage("CoolPropTPP", None, "")
        ctpp.LoadCoolProp(["Water", ],  # compound IDs in the mixture
                          [],  # unused (compund CAS numbers)
                          {"Liquid": eLiquid},  # dictionary {'phaseLabel' : stateOfAggregation}
                          eMass,  # optional: default basis (eMole by default)
                          {"backend": "HEOS",
                           "referenceState": "DEF",
                           "debugLevel": "0"})  # optional: options dictionary (empty by default)


        # Equations
        eq = self.CreateEquation("rho", "Density calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.rho(x) - self.tpp.rho(self.PR()*self.P(x),self.T(),[1.0,])

        eq = self.CreateEquation("mu", "Viscosity calculation.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.mu(x) - self.tpp.mu(self.PR()*self.P(x),self.T(),[1.0,])

        eq = self.CreateEquation("kappa", "Thermal conductivity calculation.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.kappa(x) - self.tpp.kappa(self.PR()*self.P(x),self.T(),[1.0,])

        eq = self.CreateEquation("v", "Velocity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D() ** 2
        eq.Residual = self.v(x) - self.k() / ( self.rho(x) * A )

        eq = self.CreateEquation("f_D", "Friction Factor")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Re = self.D() * self.v(x) * self.rho(x) / self.mu(x)
        eq.Residual = 1 / Sqrt(self.fD(x)) + 2. * Log10(self.ep() / 3.7 / self.D() + 2.51 / Re / Sqrt(self.fD(x)))

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        hL = 0.5 * self.fD(x) * self.v(x) ** 2 / ( self.D() * self.g() )
        DeltaP = self.g() * self.rho(x) * self.L() * hL / self.PR()
        eq.Residual = dP_dx(x) + DeltaP




class sim_test1b(daeSimulation):

    def __init__(self):

        daeSimulation.__init__(self)

        self.m = PipeWithTPP("test1b")
        self.m.Description = "Testing the solution for a pipe with Thermal Package"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):


        self.m.x.CreateStructuredGrid(10, 0.0, 1.0)

        # Setting Parameter values
        self.m.T.SetValue( 300 * K)
        self.m.g.SetValue( 9.81 * m/ s**2 )
        self.m.D.SetValue( 4.026*0.0254 * m )
        self.m.L.SetValue( 1.0 * m )
        #self.m.rho.SetValue( 1000 * kg / m**3 )
        #self.m.mu.SetValue( 0.001 * Pa * s )
        self.m.ep.SetValue( 0.0018*0.0254 * m )
        self.m.PR.SetValue( 100000. * Pa )


    def SetUpVariables(self):

        P0 = 1.0
        P1 = 0.995

        # Setting Variable Initial Guesses
        self.m.fD.SetInitialGuesses(0.018 * unit())
        self.m.v.SetInitialGuesses(2. * m / s)
        self.m.rho.SetInitialGuesses(1000. * kg / m**3 )
        self.m.mu.SetInitialGuesses(0.001 * Pa * s )
        self.m.k.SetInitialGuesses(10. * kg / s)
        self.m.P.AssignValue(0, P0 * unit())
        self.m.P.AssignValue(10, P1 * unit())

        Nx = self.m.x.NumberOfPoints

        for i in range(Nx):
            self.m.P.SetInitialGuess(i, P0 + i / 10 * (P0-P1) * unit())

        pass


    def Run(self):
        cfg = daeGetConfig()
        cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)



def test_test1a():

    data = main(simulation = sim_test1a())

    assert round( data['v']['Values'][0] - 2.36909, 4) == 0.0

def test_test1b():

    data = main(simulation = sim_test1b())

    assert round( data['v']['Values'][0][0] - 2.36909, 4) == 0.0


def test_test2():

    data = main(simulation = sim_test2())

    print(data)

    assert round( data['pipe_0.v']['Values'][0] - 2.36909, 4) == 0.0


def test_test3():

    data = main(simulation = sim_test3())

    print(data)

    assert round( data['pipe_0.v']['Values'][0] - 1.64589, 4) == 0.0



class sim_test4(daeSimulation):


    def __init__(self):

        import models.pipe as models_pipe
        daeSimulation.__init__(self)

        self.m = models_pipe.Pipe("test4")
        self.m.Description = "Testing the solution for a pipe"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        self.m.x.CreateStructuredGrid(10, 0.0, 1.0)
        # Setting Parameter values
        self.m.D.SetValue( 4.026*0.0254 * m )
        self.m.L.SetValue( 1.0 * m )
        self.m.ep.SetValue( 0.0018*0.0254 * m )

    def SetUpVariables(self):

        P0 = 1.0
        P1 = 0.995

        # Setting Variable Initial Guesses
        self.m.fD.SetInitialGuesses(0.03 * unit())
        self.m.v.SetInitialGuesses(2. * m / s)
        self.m.rho.SetInitialGuesses(1000. * kg / m**3 )
        self.m.mu.SetInitialGuesses(0.001 * Pa * s )
        # self.m.k.SetInitialGuesses(10. * kg / s)
        self.m.P.AssignValue(0, P0 * unit())
        self.m.P.AssignValue(10, P1 * unit())
        self.m.T.AssignValue(0, 300. * K)
        self.m.k.SetInitialGuess(0, 10. * kg / s)

        Nx = self.m.x.NumberOfPoints

        for i in range(Nx):
            self.m.P.SetInitialGuess(i, P0 + i / 10 * (P0-P1) * unit())

    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)



def test_test4():

    data = main(simulation = sim_test4())

    print(data)

    assert round( data['v']['Values'][0] - 1.64589, 4) == 0.0



def main(simulation = None):

    #log = daeStdOutLog()
    log = daePythonStdOutLog()
    solver = daeIDAS()

    cfg = daeGetConfig()
    cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')
    cfg.SetString('daetools.core.printInfo', 'true')
    cfg.SetString('daetools.IDAS.printInfo', 'true')

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

    methods = {'test1a': sim_test1a, 'test1b': sim_test1b, 'test2': sim_test2, 'test3': sim_test3, 'test4': sim_test4, }

    parser = argparse.ArgumentParser(description='Test pipe.')
    parser.add_argument('case',type=str,
                        help='Simulation case')

    args = parser.parse_args()

    data = main(simulation = methods[args.case]())

    pprint.pprint(data, indent = 4)
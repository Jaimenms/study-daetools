
# !/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import inspect
import pprint
import argparse

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import daeDelegateDataReporter, daeJSONFileDataReporter


# Standard variable types are defined in variable_types.py
from pyUnits import m, g, kg, s, K, Pa, mol, J, W, rad


class Nusselt(daeModel):

    def __init__(self, Name, Parent=None, Description=""):

        daeModel.__init__(self, Name, Parent, Description)

        no_pos_t = daeVariableType("no_pos_t", dimless, 1e-10, 1e10, 1.0, 1e-06)
        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-3, 1e+02, 0.1, 1e-05)
        specific_heat_capacity_t = daeVariableType("specific_heat_capacity_t", (J ** (1))*(K ** (-1))*(kg ** (-1)), 1e-10, 1e+10, 1000., 1e-05)

        self.D = daeParameter("D", m, self, "External pipe diameter")

        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Fluid thermal conductivity")
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid")
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid")
        self.cp = daeVariable("cp", specific_heat_capacity_t, self, "Cp")

        self.v = daeVariable("v", velocity_t, self, "Fluid Velocity")

        self.Re = daeVariable("Re", no_t, self, "Reynolds")

        self.prandtl = daeVariable("prandtl", no_pos_t, self, "Prandtl")
        self.nusselt = daeVariable("nusselt", no_pos_t, self, "Nusselt")


    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("Prandtl", "Prandtl")
        eq.Residual = self.prandtl() - self.cp() * self.mu() / self.kappa()

        eq = self.CreateEquation("Re", "Re")
        eq.Residual = self.Re() - self.D() * self.v() * self.rho() / self.mu()

        eq = self.CreateEquation("Nusselt", "Nusselt")
        eq.Residual = self.nusselt() - 0.027 * Pow(self.Re(), 4/5) * Pow(self.prandtl(), 1/3)



class sim_01(daeSimulation):


    def __init__(self):

        import models.pipe as models_pipe
        daeSimulation.__init__(self)

        self.m = Nusselt("test")
        self.m.Description = "Testing"

        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):

        self.m.D.SetValue( 4.026*0.0254 * m )

    def SetUpVariables(self):

        self.m.kappa.AssignValue(0.6094991151038377 * W/(K * m))
        self.m.mu.AssignValue(0.0008537426062537152 * Pa * s)
        self.m.rho.AssignValue(996.5563403888981 * kg / m ** 3)
        self.m.cp.AssignValue(4180.639522022376 * J/(K * kg))

        self.m.Re.SetInitialGuess(100000. * unit() )
        self.m.prandtl.SetInitialGuess(5. * unit() )
        self.m.v.AssignValue(2.39 * m / s)
        #self.m.nusselt.AssignValue(1125.779997967431 * unit() )
        self.m.nusselt.SetInitialGuess(1000. * unit() )

    def Run(self):
        # A custom operating procedure, if needed.
        # Here we use the default one:
        daeSimulation.Run(self)



def test_case():

    data = main(simulation = sim_01())

    print(data)

    assert round( data['v']['Values'][0][0] - 2.38928, 4) == 0.0


def main(simulation = sim_01()):

    # cfg = daeGetConfig()
    # cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')
    #cfg.SetString('daetools.core.printInfo', 'true')
    #cfg.SetString('daetools.IDAS.printInfo', 'true')

    log = daeStdOutLog()
    solver = daeIDAS()

    simulation.m.SetReportingOn(True)

    datareporter = daeDelegateDataReporter()

    dr = daeJSONFileDataReporter()

    datareporter.AddDataReporter(dr)

    json_filename = simulation.report_filename
    simName = os.path.basename(json_filename) + '_' + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    dr.Connect(json_filename, simName)

    print('%s (%s): %s' % (
    dr.__class__.__name__, dr.ConnectString, 'connected' if dr.IsConnected() else 'NOT connected'))

    simulation.Initialize(solver, datareporter, log)
    simulation.SolveInitial()
    simulation.Finalize()

    return json.load(open(simulation.report_filename))


if __name__ == "__main__":

    data = main()

    pprint.pprint(data, indent = 4)
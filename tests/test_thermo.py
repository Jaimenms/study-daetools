
# !/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import inspect
import pprint

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from daetools.pyDAE import *
from daetools_reporters import setupDataReporters


# Standard variable types are defined in variable_types.py
from pyUnits import m, g, kg, s, K, Pa, mol, J, W

def test_coolprop_single_phase():
    """
    Test the CoolProp library.
    Default backend is HEOS and the default reference state is DEF.
    Different backends/ref.states can be specified in the options dictionary.
    Nota bene:
      CoolProp does not have transport properties models defined for some compounds (i.e. CarbonMonoxide).
    """
    ctpp = daeThermoPackage("CoolPropTPP", None, "")
    ctpp.LoadCoolProp(["Water",],                      # compound IDs in the mixture
                      [],                              # unused (compund CAS numbers)
                      {"Liquid"  : eLiquid},           # dictionary {'phaseLabel' : stateOfAggregation}
                      eMass,                           # optional: default basis (eMole by default)
                      {"backend"        : "HEOS",
                       "referenceState" : "DEF",
                       "debugLevel"     : "0"})        # optional: options dictionary (empty by default)


    P = 1.0e5 # Must be in Pascals
    T = 300   # Must be in Kelvins
    x = [1,]

    assert round(996.5563 - ctpp._CalcSinglePhaseScalarProperty("density",  P, T, x, "Liquid", eMass), 3) == 0.
    assert round(4180.6395 - ctpp._CalcSinglePhaseScalarProperty("heatCapacityCp",  P, T, x, "Liquid", eMass), 3) == 0.
    assert round(0.6095 - ctpp._CalcSinglePhaseScalarProperty("thermalConductivity",  P, T, x, "Liquid", eUndefinedBasis), 3) == 0.
    assert round(112653.6796 - ctpp._CalcSinglePhaseScalarProperty("enthalpy",  P, T, x, "Liquid", eMass), 3) == 0.
    assert round(0.000853742 - ctpp._CalcSinglePhaseScalarProperty("viscosity",  P, T, x, "Liquid", eUndefinedBasis), 3) == 0.


class modTutorial(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self, Name, Parent, Description)

        self.tpp = daeThermoPackage("TPP", self, "")
        # CapeOpen thermo packages are available only in Windows.
        # For other platforms use CoolProp thermo package.
        self.tpp.LoadCoolProp(["Water",],       # compound IDs in the mixture
            [],                                 # compund CAS numbers (optional)
            {'Liquid': eLiquid},                # dictionary {'phaseLabel' : stateOfAggregation}
            eMole,                              # default basis is eMole (other options are eMass or eUndefinedBasis)
            {})                                 # options dictionary (defaut is empty)

        self.m     = daeParameter("m", kg, self, "Mass of water")

        self.cp    = daeVariable("c_p",   specific_heat_capacity_t, self, "Specific heat capacity of the liquid")
        self.kappa = daeVariable("k",     thermal_conductivity_t,   self, "Thermal conductivity of the liquid")
        self.mu    = daeVariable("mu",    dynamic_viscosity_t,      self, "Viscosity of the liquid")
        self.Qin   = daeVariable("Q_in",  power_t,                  self, "Power of the heater")
        self.T     = daeVariable("T",     temperature_t,            self, "Temperature of the liquid")

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # P, T, x arguments can be floats or adouble objects/arrays as illustrated here.
        P = 1e5 * Pa
        T = self.T()
        x = [1.0,]

        # daeThermoPackage uses the default basis to calculate the properties.
        # If a property for a different basis is required it should be specified as a keyword argument 'basis'=eMole|eMass|eUndefinedBasis.
        # Also, it needs the phase information. If there is only one phase it will be used as a default.
        # Otherwise, it must be specified as a keyword argument 'phase'='phaseLabel'.
        # For instance:
        #    self.tpp.kappa(P,T,x, phase='Vapor', basis=eUndefinedBasis)

        # Calculate transport properties:
        # 1. Specific heat capacity (cp) in J/(kg*K)
        #    TPP package returns cp in J/(mol.K) (mole basis) - for some reasons it refuses to calculate it for the mass basis.
        #    Therefore, divide it by the molar mass of the mixture.
        eq = self.CreateEquation("C_p", "Specific heat capacity as a function of the temperature.")
        eq.Residual = self.cp() - self.tpp.cp(P ,T ,x) / self.tpp.M(P ,T ,x)

        # 2. Thermal conductivity (kappa) in W/(m*K)
        eq = self.CreateEquation("kappa", "Thermal conductivity as a function of the temperature.")
        eq.Residual = self.kappa() - self.tpp.kappa(P ,T ,x)

        # 3. Dynamic viscosity (mu) in Pa*s
        eq = self.CreateEquation("mu", "Viscosity as a function of the temperature.")
        eq.Residual = self.mu() - self.tpp.mu(P ,T ,x)

        # Simple integral heat balance for the liquid
        eq = self.CreateEquation("HeatBalance", "Integral heat balance equation")
        eq.Residual = self.m() * self.cp() * dt(self.T()) - self.Qin()


class sim_test1(daeSimulation):

    def __init__(self):
        daeSimulation.__init__(self)
        self.m = modTutorial("test")
        self.m.Description = ""
        self.report_filename = __file__ + '.json'


    def SetUpParametersAndDomains(self):
        self.m.m.SetValue(1 * kg)

    def SetUpVariables(self):
        self.m.Qin.AssignValue(500 * W)
        self.m.T.SetInitialCondition(283 * K)

def run(**kwargs):
    # Thermo-physical property packages are not supported by the Compute Stack approach.
    # Therefore, activate the old approach.
    cfg = daeGetConfig()
    cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')

    simulation = sim_test1()
    return daeActivity.simulate(simulation, reportingInterval = 10,
                                timeHorizon       = 500,
                                **kwargs)


def test_test1():

    data = main(simulation = sim_test1())

    assert round( data['c_p']['Values'][0] - 4195.4086, 3) == 0.0


def main(simulation = sim_test1()):

    cfg = daeGetConfig()

    cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')

    log = daeStdOutLog()
    solver = daeIDAS()

    simulation.m.SetReportingOn(True)

    datareporter = setupDataReporters(simulation)

    simulation.Initialize(solver, datareporter, log)

    # Solve at time = 0 (initialization)
    simulation.SolveInitial()

    # Clean up
    simulation.Finalize()

    return json.load(open(simulation.report_filename))


if __name__ == "__main__":

    pprint.pprint(main(), indent = 4)
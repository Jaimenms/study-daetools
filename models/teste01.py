#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
***********************************************************************************
                            tutorial1.py
                DAE Tools: pyDAE module, www.daetools.com
                Copyright (C) Dragan Nikolic
***********************************************************************************
DAE Tools is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License version 3 as published by the Free Software
Foundation. DAE Tools is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with the
DAE Tools software; if not, see <http://www.gnu.org/licenses/>.
************************************************************************************
"""
__doc__ = """
A ser preenchido
"""

import os, sys, tempfile
from time import localtime, strftime
from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from fluids.friction import friction_factor as ff

# Standard variable types are defined in variable_types.py
from pyUnits import m, kg, s, K, Pa, mol, J, W


class extfn_DarcyFactor(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, Re, eD):
        arguments = {}
        arguments["Re"] = Re
        arguments["eD"] = eD

        # Instantiate interp1d object and initialize interpolation using supplied (x,y) values
        self.ff_func = ff

        # During the solver iterations, the function is called very often with the same arguments
        # Therefore, cache the last interpolated value to speed up a simulation
        self.cache = None

        # Counters for performance (just an info; not really needed)
        self.counter = 0
        self.cache_counter = 0

        daeScalarExternalFunction.__init__(self, Name, Model, units, arguments)

    def Calculate(self, values):
        # Increase the call counter every time the function is called
        self.counter += 1

        # Get the argument from the dictionary of arguments' values.
        Re = values["Re"].Value
        eD = values["eD"]

        # First check if an interpolated value was already calculated during the previous call
        # If it was return the cached value (derivative part is always equal to zero in this case)
        if self.cache:
            if self.cache[0] == Re and self.cache[1] == eD:
                self.cache_counter += 1
                return adouble(self.cache[2])

        # The time received is not in the cache and has to be interpolated.
        # Convert the result to float datatype since daetools can't accept
        # numpy.float64 types as arguments at the moment.
        res = self.ff_func(Re=Re, eD=eD)

        # Save it in the cache for later use
        self.cache = (Re, eD, ff)

        return res


class modTutorial(daeModel):
    def __init__(self, Name, Parent=None, Description=""):
        daeModel.__init__(self, Name, Parent, Description)

        # Distribution domain is a general term used to define an array of different objects (parameters, variables, etc).
        # daeDomain constructor accepts three arguments:
        #  - Name: string
        #  - Parent: daeModel object (indicating the model where the domain will be added)
        #  - Units: unit object (specifying the units, such as 'm' for we are modelling a space distribution)
        #  - Description: string (optional argument; the default value is an empty string)
        # All naming conventions (introduced in Whats_the_time example) apply here as well.
        # Again, domains have to be declared as members of the new model class (like all the other objects)
        self.x = daeDomain("x", self, m / m, "X axis domain")

        # Parameter can be defined as a time invariant quantity that will not change during a simulation.
        # daeParameter constructor accepts three arguments:
        #  - Name: string
        #  - Units: unit object (specifying the units, such as 'kg', 'W', 'm/s', ...)
        #  - Parent: daeModel object (indicating the model where the domain will be added)
        #  - Description: string (optional argument; the default value is an empty string)
        # All naming conventions (introduced in whats_the_time example) apply here as well.
        self.L = daeParameter("L_tube", m, self, "Tube length")
        self.D = daeParameter("d_tube", m, self, "Tube internal diameter")
        self.eD = daeParameter("&epsilon;_D", m / m, self, "Tube roughness / diameter")
        self.A = daeParameter("A_tube", m ** 2, self, "Liquid cross section area")
        self.q0 = daeParameter("q0", kg / s, self, "Referential mass flow rate")
        self.Ta = daeParameter("Ta", K, self, "Temperature at inlet")
        self.Pa = daeParameter("Pa", Pa, self, "Pressure at inlet")
        self.Pb = daeParameter("Pb", Pa, self, "Pressure at outlet")

        self.rhoL = daeVariable("&rho;_L", density_t, self)
        self.rhoL.DistributeOnDomain(self.x)
        self.rhoL.Description = "Liquid Density in kg/m3"

        self.cpL = daeVariable("c_pL",specific_heat_capacity_t, self)
        self.cpL.DistributeOnDomain(self.x)
        self.cpL.Description = "Heat capacity in J/(kg*K)"

        self.alphaL = daeVariable("&alpha;_L", fraction_t, self)
        self.alphaL.DistributeOnDomain(self.x)
        self.alphaL.Description = "Liquid hold-up"

        self.v = daeVariable("v", velocity_t, self)
        self.v.DistributeOnDomain(self.x)
        self.v.Description = "Velocity of the liquid"

        self.qL = daeVariable("q_L", no_t, self)
        self.qL.DistributeOnDomain(self.x)
        self.qL.Description = "Mass flow rate of the liquid (kg/s)"

        self.T = daeVariable("T", temperature_t, self)
        self.T.DistributeOnDomain(self.x)
        self.T.Description = "Temperature of the liquid"

        self.P = daeVariable("P", pressure_t, self)
        self.P.DistributeOnDomain(self.x)
        self.P.Description = "Pressure of the liquid"

        self.Re = daeVariable("Re", no_t, self)
        self.Re.DistributeOnDomain(self.x)
        self.Re.Description = "Reynolds of the liquid"

        self.ff = daeVariable("f", no_t, self)
        self.ff.DistributeOnDomain(self.x)
        self.ff.Description = "Friction Factor"

        self.tpp = daeThermoPackage("TPP", self, "")
        self.tpp.LoadCoolProp(["Water",],  # compound IDs in the mixture
                              [],  # compund CAS numbers (optional)
                              {'Liquid': eLiquid},  # dictionary {'phaseLabel' : stateOfAggregation}
                              eMass,  # default basis is eMole (other options are eMass or eUndefinedBasis)
                              {})  # options dictionary (defaut is empty)

    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        rhoL = lambda x: self.rhoL(x)
        cpL = lambda x: self.cpL(x)
        alphaL = lambda x: self.alphaL(x)
        v = lambda x: self.v(x)
        qL = lambda x: self.qL(x)
        T = lambda x: self.T(x)
        P = lambda x: self.P(x)
        Re = lambda x: self.Re(x)
        ff = lambda x: self.ff(x)

        # Friction Factor
        #N = self.x.NumberOfPoints
        #for i in range(N):
        #    eq = self.CreateEquation("DarcyFrictionFactor%d" % i)
        #    eq.Residual = ff(i) - extfn_DarcyFactor("FrictionFactor", self, m / m, self.Re(i), self.eD())

        N = self.x.NumberOfPoints
        F_ff = numpy.zeros(N, dtype=object)
        F_ff[:] = [extfn_DarcyFactor("FrictionFactor", self, m / m, self.Re(i), self.eD()) for i in range(N)]
        for i in range(N):
            eq = self.CreateEquation("Darcy Friction Factor(%d)" % i)
            eq.Residual = ff(i) - F_ff[i]

        #eq = self.CreateEquation("ff", "Darcy Friction Factor as a function of velociy, temperature and pressure.")
        #x = eq.DistributeOnDomain(self.x, eClosedClosed)
        #eq.Residual = ff(x) - F_ff

        # Reynolds (Re)
        eq = self.CreateEquation("Re", "Reynolds as a function of velociy, temperature and pressure.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = Re(x) - self.D() * v(x) * rhoL(x) / self.tpp.mu(P(x), T(x), [1, ])

        # Density (rho) in kg/m**3
        eq = self.CreateEquation("rho_L", "Density as a function of temperature and pressure.")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = rhoL(x) - self.tpp.rho(P(x),T(x),[1,])

        # Calculate transport properties:
        # Specific heat capacity (cpL) in J/kg
        eq = self.CreateEquation("c_pL", "Heat capacity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = cpL(x) - self.tpp.cp(P(x),T(x),[1,])

        # Friction Factor

        # Mass Flow rate (q_L) in kg/s
        eq = self.CreateEquation("q_L", "Mass flow rate")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = qL(x) - v(x) * self.A() * alphaL(x) * rhoL(x) / self.q0()

        eq = self.CreateEquation("Biofilm", "Biofilm formation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        #eq.Residual = dt(alphaL(x)) - self.Salfa(x)
        eq.Residual = dt(alphaL(x))

        eq = self.CreateEquation("MassBalance", "Mass balance equation valid on open x domain")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = dt(rhoL(x)*alphaL(x)) +\
                      d(rhoL(x)*self.v(x)*alphaL(x),self.x,eCFDM)/self.L()

        eq = self.CreateEquation("MomentumBalance", "Momentum balance equation valid on open x domain")
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(rhoL(x)*v(x)*alphaL(x)) +\
                      d(rhoL(x)*v(x)**2,self.x,eCFDM)/self.L() +\
                      alphaL(x)*d(P(x),self.x,eCFDM)/self.L() +\
                      ff(x)/self.D()*self.rhoL()*v(x)**2/2

        eq = self.CreateEquation("HeatBalance", "Heat balance equation valid on open x domain")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)
        eq.Residual = dt(rhoL(x)*cpL(x)*T(x)*alphaL(x)) +\
                      d(rhoL(x)*cpL(x)*T(x)*v(x),self.x,eCFDM)/self.L()# - self.Sheat(x)

        # Boundary conditions are treated as ordinary equations, and the special eLowerBound and eUpperBound flags
        # are used to define the position of the boundary.
        eq = self.CreateEquation("BC_inlet_T", "Dirichlet boundary condition at the tube inlet")
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = T(x) - self.Ta()

        eq = self.CreateEquation("BC_inlet_P", "Dirichlet boundary condition at the tube inlet")
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = P(x) - self.Pa()

        eq = self.CreateEquation("BC_outlet_P", "Dirichlet boundary condition at the tube outlet")
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = P(x) - self.Pb()



class simTutorial(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = modTutorial("teste01")
        self.m.Description = __doc__

    def SetUpParametersAndDomains(self):
        # The function CreateStructuredGrid can be used to create a structured grid. It accepts 3 arguments:
        #  - NoIntervals: integer
        #  - LowerBound:  float
        #  - UpperBound:  float
        # Here we create 2D structured gird with 10x10 intervals in x and y domains
        # (resulting in 11x11 points in both directions).
        self.m.x.CreateStructuredGrid(5, 0, 1)

        self.m.L.SetValue(5 * m)
        self.m.D.SetValue(3 / 4 * 0.0254 * m)
        self.m.eD.SetValue(0.00004)
        self.m.A.SetValue(0.25*3.14*(3/4*0.0254)**2 * (m ** 2))
        self.m.Ta.SetValue(300 * K)
        self.m.q0.SetValue(1 * kg / s)
        self.m.Pa.SetValue(4e5 * Pa)
        self.m.Pb.SetValue(3.8e5 * Pa)

        self.InitialConditionMode = eQuasiSteadyState


def setupDataReporters(simulation):

    datareporter = daeDelegateDataReporter()

    dr = daeJSONFileDataReporter()

    datareporter.AddDataReporter(dr)

    # Connect data reporters
    modelName = simulation.m.Name
    simName = modelName + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())

    #directory = tempfile.gettempdir()
    directory = "./"

    json_filename = os.path.join(directory, "%s.json" % modelName)
    dr.Connect(json_filename, simName)

    print('%s (%s): %s' % (
    dr.__class__.__name__, dr.ConnectString, 'connected' if dr.IsConnected() else 'NOT connected'))

    return datareporter


def run(**kwargs):
    simulation = simTutorial()
    return daeActivity.simulate(simulation, reportingInterval=5,
                                timeHorizon=500,
                                datareporter=setupDataReporters(simulation),
                                **kwargs)


if __name__ == "__main__":
    guiRun = False if (len(sys.argv) > 1 and sys.argv[1] == 'console') else True
    run(guiRun=guiRun)


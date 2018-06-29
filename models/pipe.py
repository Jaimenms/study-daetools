__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
# These are methods and constants imported from daeTools C++ code
daeVariable = daeVariable
daeVariableType = daeVariableType
daeDomain = daeDomain
Constant = Constant
daeParameter = daeParameter
eClosedClosed = eClosedClosed
eOpenClosed = eOpenClosed
unit = unit
Abs = Abs
Log = Log
Log10 = Log10
Exp = Exp
Sin = Sin
eUpperBound = eUpperBound
eLowerBound = eLowerBound
eCFDM = eCFDM
dt = dt
d = d
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.edge import Edge
except:
    from .edge import Edge

from water_properties import density, viscosity, conductivity, heat_capacity


class Pipe(Edge):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        Edge.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_domains(self):

        Edge.define_domains(self)

        self.x = daeDomain("x", self, unit(), "X axis domain")
        self.domain_len = self.data['domains']['N']


    def define_constants(self):

        Edge.define_constants(self)

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa) # bar


    def define_parameters(self):

        Edge.define_parameters(self)

        self.tetha = daeParameter("tetha", rad, self, "Angle")
        self.Di = daeParameter("Di", m, self, "Inside pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.Kub = daeParameter("Kub", unit(), self, "Upper bound concentrated accident coefficient")
        self.Klb = daeParameter("Klb", unit(), self, "Lower bound concentrated accident coefficient")
        self.Npipes = daeParameter("Npipes", unit(), self, "Number of pipes")


    def define_variables(self):

        Edge.define_variables(self)

        # Variable types
        diameter_t = daeVariableType("diameter_t", (m ** (1)), 0.00001, 10.0, 0.5, 1e-05)
        darcy_t = daeVariableType("darcy_t", dimless, 0.0001, 5, 0.018, 1e-03)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-6, 100.0, 0.1, 1e-06)
        deltap_t = daeVariableType("deltap_t", Pa / m, -1e6, 1e6, 2000, 1e-05)

        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-6, 10.0, 1.0, 1e-06)

        # Secondary variables
        self.H = daeVariable("H", power_t, self, "Fluid Entalphy", [self.x, ])
        self.fD = daeVariable("fD", darcy_t, self, "Darcy friction factor", [self.x, ])
        self.fDturb = daeVariable("fDturb", darcy_t, self, "Turbulent Darcy friction factor", [self.x, ])
        self.fDlam = daeVariable("fDlam", darcy_t, self, "Laminar Darcy friction factor", [self.x, ])
        self.fregime = daeVariable("fregime", no_t, self, "Flow regime selector", [self.x, ])
        self.D = daeVariable("D", diameter_t, self, "Internal flow diameter", [self.x, ])
        self.v = daeVariable("v", velocity_t, self, "Internal flow velocity", [self.x, ])
        self.dp = daeVariable("dp", deltap_t, self, "Delta Pressure", [self.x, ])

        # Fluid Properties
        self.Re = daeVariable("Re", no_t, self, "Viscosity of the liquid", [self.x, ])
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid", [self.x, ])
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid", [self.x, ])
        self.cp = daeVariable("cp", specific_heat_capacity_t, self, "Heat capacity of the liquid", [self.x, ])
        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Thermal Conductivity of the liquid", [self.x, ])

        # State variables
        self.Pub = daeVariable("Pub", pressure_t, self, "Upper Bound Pressure")
        self.Plb = daeVariable("Plb", pressure_t, self, "Lower Bound Pressure")
        self.P = daeVariable("P", pressure_t, self, "Fluid Pressure", [self.x, ])
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", [self.x, ])
        self.Tub = daeVariable("Tub", water_temperature_t, self, "Upper Bound Fluid Temperature")
        self.Tlb = daeVariable("Tlb", water_temperature_t, self, "Lower Bound Fluid Temperature")

        # Concentrated State Variables
        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate")


    def eq_velocity(self):

        eq = self.CreateEquation("v", "v_fluid_velocity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        A = 0.25 * 3.14 * self.D(x) ** 2
        eq.Residual = self.v(x) - self.k() / ( self.rho(x) * A )


    def eq_entalphy(self):

        eq = self.CreateEquation("Entalphy", "Fluid Entalphy")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.k() * self.cp(x) * self.T(x) - self.H(x)


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.D(x) - self.Di()


    def eq_Reynolds(self):

        eq = self.CreateEquation("Re", "Reynolds")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.Re(x) - self.D(x) * Abs(self.v(x)) * self.rho(x) / self.mu(x)


    def eq_fregime(self):

        eq = self.CreateEquation("fregime", "Flow regime selector")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Re = self.Re(x)
        fregime = (1 + Exp(-10 * (Re / 3150 - 1))) ** -1
        eq.Residual = self.fregime(x) - fregime


    def eq_fDlaminar(self):

        eq = self.CreateEquation("fDlam", "Laminar darcy factor")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Re = self.Re(x)
        fDlam = 64/Re
        eq.Residual = self.fDlam(x) - fDlam


    def eq_fDturb(self):

        eq = self.CreateEquation("fDturb", "Turbulent darcy factor")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Re = self.Re(x)
        ep = self.ep()
        D = self.D(x)
        A = (2.457 * Log(((7 / Re) ** 0.9 + 0.27 * ep / D) ** -1)) ** 16
        B = (37530 / Re) ** 16
        ff = 2 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12)
        fDturb  = 4 * ff
        eq.Residual = self.fDturb(x) - fDturb


    def eq_fD(self):

        self.stnFlowRegimeSetup = self.STN("FlowRegimeSetup")

        self.STATE("Laminar")

        eq = self.CreateEquation("fDForcedLaminar", "Darcy factor - Laminar")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        fDlam = self.fDlam(x)
        eq.Residual = self.fD(x) - fDlam

        # self.ON_CONDITION(Time() > Constant(100000000*s), switchToStates     = [ ('FlowRegimeSetup', 'Calculated') ],
        #                                             setVariableValues  = [],
        #                                             triggerEvents      = [],
        #                                             userDefinedActions = [] )


        self.STATE("Turbulent")

        eq = self.CreateEquation("fDForcedTurbulent", "Darcy factor - Turbulent")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        fDturb = self.fDturb(x)
        eq.Residual = self.fD(x) - fDturb

        # self.ON_CONDITION(Time() > Constant(100*s), switchToStates     = [ ('FlowRegimeSetup', 'Calculated') ],
        #                                             setVariableValues  = [],
        #                                             triggerEvents      = [],
        #                                             userDefinedActions = [] )


        self.STATE("Calculated")

        eq = self.CreateEquation("fDCalculated", "Darcy factor - Calculated")

        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        fDlam = self.fDlam(x)
        fDturb = self.fDturb(x)
        fregime = self.fregime(x)
        fD = self.fD(x) - ((fDturb - fDlam) * fregime + fDlam)
        eq.Residual = self.fD(x) - fD

        self.END_STN()


    def eq_delta_pressure(self):

        eq = self.CreateEquation("DeltaPressure", "Delta Pressure")
        x = eq.DistributeOnDomain(self.x,eClosedClosed)
        Re = self.Re(x)
        fD = self.fD(x) #TODO
        hL = 0.5 * fD * Abs(self.v(x)) * self.v(x) / ( self.D(x) * self.g )
        DeltaP = self.g * self.rho(x) * hL
        eq.Residual = self.dp(x) - DeltaP


    def eq_lowerbound_pressure(self):

        eq = self.CreateEquation("LowerBoundPressure", "Lower Bound Pressure")
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        DeltaPlb = 0.5 * self.Klb() * self.rho(x) * self.v(x) ** 2
        eq.Residual = self.P(x) + DeltaPlb - self.Plb()


    def eq_upperbound_pressure(self):

        eq = self.CreateEquation("UpperBoundPressure", "Upper Bound Pressure")
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        DeltaPlub = 0.5 * self.Kub() * self.rho(x) * self.v(x) ** 2
        eq.Residual = self.P(x) - DeltaPlub - self.Pub()


    def eq_lowerbound_temperature(self):

        eq = self.CreateEquation("LowerBoundPressure", "Lower Bound Pressure")
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.T(x) - self.Tlb()


    def eq_upperbound_temperature(self):

        eq = self.CreateEquation("UpperBoundPressure", "Upper Bound Pressure")
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = self.T(x) - self.Tub()


    def eq_mommentum_balance(self):

        eq = self.CreateEquation("MomBal", "Momentum balance")
        x = eq.DistributeOnDomain(self.x,eOpenClosed)

        A = 0.25 * 3.14 * self.D(x) ** 2

        ## Nessa forma é Darcy TODO Corrigir isso!!!
        hL = 0.5 * self.fD(x) * Abs(self.v(x)) * self.v(x) / ( self.D(x) * self.g )
        DeltaP = self.g * self.rho(x) * hL

        #eq.Residual = dt(A * self.rho(x) * self.v(x)) + self.k() * d(self.v(x), self.x, eCFDM) / self.L() + d(
        #    A * self.P(x), self.x, eCFDM) / self.L() + A * DeltaP + A * self.rho(x) * self.g * Sin(self.tetha())
        eq.Residual =  dt(A * self.rho(x) * self.v(x)) + self.k() * d(self.v(x), self.x, eCFDM) / self.L() + A * d(self.P(x), self.x, eCFDM) / self.L()  +  A * DeltaP +  A * self.rho(x) * self.g * Sin(self.tetha())

    def eq_heat_balance(self):

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        x = eq.DistributeOnDomain(self.x, eOpenClosed)

        A = 0.25 * 3.14 * self.D(x) ** 2
        eq.Residual = self.rho(x) * self.cp(x) * dt(A * self.T(x)) + self.k() * self.cp(x) * d( self.T(x), self.x, eCFDM) / self.L()


    def eq_water_properties(self):

        self.stnWaterPropertiesSetup = self.STN("WaterPropertiesSetup")

        self.STATE("Fixed")

        eq = self.CreateEquation("Density", "Density calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        rho = density(300., 100000., simplified=True)
        eq.Residual = self.rho(x) - rho * Constant(1 * (kg ** (1))*(m ** (-3)))

        eq = self.CreateEquation("Viscosity", "Viscosity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        mu = viscosity(300., 100000., simplified=True)
        eq.Residual = self.mu(x) - mu * Constant(1 * (Pa ** (1))*(s ** (1)))

        eq = self.CreateEquation("ThermalConductivity", "Thermal Conductivity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        kappa = conductivity(300., 100000., simplified=True)
        eq.Residual = self.kappa(x) - kappa * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))

        eq = self.CreateEquation("HeatCapacity", "Heat Capacity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        cp = heat_capacity(300., 100000., simplified=True)
        eq.Residual = self.cp(x) - cp * Constant(1 * (J ** (1))*(K ** (-1))*(kg ** (-1)))

        self.STATE("Variable")

        eq = self.CreateEquation("Density", "Density calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        rho = density( self.T(x) / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified = True)
        eq.Residual = self.rho(x) - rho * Constant(1 * (kg ** (1))*(m ** (-3)))

        eq = self.CreateEquation("Viscosity", "Viscosity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        mu = viscosity( self.T(x) / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified = True)
        eq.Residual = self.mu(x) - mu * Constant(1 * (Pa ** (1))*(s ** (1)))

        eq = self.CreateEquation("ThermalConductivity", "Thermal Conductivity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        kappa = conductivity( self.T(x) / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified = True)
        eq.Residual = self.kappa(x) - kappa * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))

        eq = self.CreateEquation("HeatCapacity", "Heat Capacity calculation")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        cp = heat_capacity( self.T(x) / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified = True)
        eq.Residual = self.cp(x) - cp * Constant(1 * (J ** (1))*(K ** (-1))*(kg ** (-1)))

        self.END_STN()


    def DeclareEquations(self):

        Edge.DeclareEquations(self)

        self.eq_fDlaminar()
        self.eq_fDturb()
        self.eq_fregime()
        self.eq_fD()
        self.eq_pressure_boundaries()
        self.eq_temperature_boundaries()
        self.eq_velocity()
        self.eq_internal_diameter()
        self.eq_entalphy()
        self.eq_mommentum_balance()
        self.eq_heat_balance()
        self.eq_Reynolds()
        self.eq_water_properties()
        self.eq_delta_pressure()
        self.eq_upperbound_pressure()
        self.eq_lowerbound_pressure()
        self.eq_upperbound_temperature()
        self.eq_lowerbound_temperature()

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

from daetools_extended.tools import daeVariable_wrapper, distribute_on_domains



class Pipe(Edge):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        Edge.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_domains(self):

        Edge.define_domains(self)

        self.x = daeDomain("x", self, unit(), "X axis domain")


    def define_constants(self):

        Edge.define_constants(self)

        self.pi = Constant( 3.14 )
        self.g = Constant(9.81 * m / s ** 2)
        self.PR = Constant(1e5 * Pa) # bar

        y_domain = []
        if hasattr(self,'y') and type(self.y) == pyCore.daeDomain:
            y_domain.append(self.y)
        self.YDomains = y_domain

        x_domain = []
        if hasattr(self,'x') and type(self.x) == pyCore.daeDomain:
            x_domain.append(self.x)
        self.XDomains = x_domain


    def define_parameters(self):

        Edge.define_parameters(self)

        self.tetha = daeParameter("tetha", rad, self, "Angle")
        self.Di = daeParameter("Di", m, self, "Inside pipe diameter")
        self.L = daeParameter("L", m, self, "Length")
        self.ep = daeParameter("epsilon", m, self, "Roughness")
        self.Kub = daeParameter("Kub", unit(), self, "Upper bound concentrated accident coefficient")
        self.Klb = daeParameter("Klb", unit(), self, "Lower bound concentrated accident coefficient")
        self.Npipes = daeParameter("Npipes", unit(), self, "Number of pipes", self.YDomains)


    def define_variables(self):

        Edge.define_variables(self)

        # Variable types
        diameter_t = daeVariableType("diameter_t", (m ** (1)), 0.00001, 10.0, 0.5, 1e-05)
        darcy_t = daeVariableType("darcy_t", dimless, 0.0001, 5, 0.018, 1e-03)
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        mass_flowrate_t = daeVariableType("mass_flowrate_t", (kg ** (1)) * (s ** (-1)), 1e-6, 100.0, 0.1, 1e-06)
        deltap_t = daeVariableType("deltap_t", Pa / m, -1e6, 1e6, 2000, 1e-05)
        velocity_t = daeVariableType("velocity_t", (m ** (1)) * (s ** (-1)), 1e-6, 10.0, 1.0, 1e-06)
        heat_per_length_t = daeVariableType("heat_per_length_t", (J ** (1)) * (m ** (-1)) * (s ** (-1)), -1e+10, 1e+10,
                                            0.1, 1e-05)

        # Secondary variables
        self.H = daeVariable("H", power_t, self, "Fluid Entalphy", self.Domains)
        self.fD = daeVariable("fD", darcy_t, self, "Darcy friction factor", self.Domains)
        self.D = daeVariable("D", diameter_t, self, "Internal flow diameter", self.Domains)
        self.v = daeVariable("v", velocity_t, self, "Internal flow velocity", self.Domains)
        self.Qout = daeVariable("Qout", heat_per_length_t, self, "Mass loss per length", self.Domains)


        # Fluid Properties
        self.Re = daeVariable("Re", no_t, self, "Viscosity of the liquid", self.Domains)
        self.mu = daeVariable("mu", dynamic_viscosity_t, self, "Viscosity of the liquid", self.Domains)
        self.rho = daeVariable("rho", density_t, self, "Density of the liquid", self.Domains)
        self.cp = daeVariable("cp", specific_heat_capacity_t, self, "Heat capacity of the liquid", self.Domains)
        self.kappa = daeVariable("kappa", thermal_conductivity_t, self, "Thermal Conductivity of the liquid", self.Domains)

        # State variables
        self.P = daeVariable("P", pressure_t, self, "Fluid Pressure", self.Domains)
        self.T = daeVariable("T", water_temperature_t, self, "Fluid Temperature", self.Domains)

        # Concentrate state variable
        self.k = daeVariable("k", mass_flowrate_t, self, "Mass flowrate", self.YDomains)

        # Bounds
        self.klb = daeVariable("klb", mass_flowrate_t, self, "Lower Bound Mass flowrate")
        self.kub = daeVariable("kub", mass_flowrate_t, self, "Upper Bound Mass flowrate")
        self.Hub = daeVariable("Hub", power_t, self, "Upper Bound Fluid Enthalpy")
        self.Hlb = daeVariable("Hlb", power_t, self, "Lower Bound Fluid Enthalpy")


    def eq_velocity(self):

        eq = self.CreateEquation("v", "v_fluid_velocity")
        xdomains = distribute_on_domains(self.XDomains, eq, eClosedClosed)
        ydomains = distribute_on_domains(self.YDomains, eq, eClosedClosed)
        domains = xdomains + ydomains

        k = daeVariable_wrapper(self.k, ydomains)

        rho = daeVariable_wrapper(self.rho, domains)
        D = daeVariable_wrapper(self.D, domains)
        v = daeVariable_wrapper(self.v, domains)

        A = 0.25 * 3.14 * D ** 2
        eq.Residual = v - k / ( rho * A )


    def eq_entalphy(self):

        eq = self.CreateEquation("Entalphy", "Fluid Entalphy")
        xdomains = distribute_on_domains(self.XDomains, eq, eClosedClosed)
        ydomains = distribute_on_domains(self.YDomains, eq, eClosedClosed)
        domains = xdomains + ydomains

        k = daeVariable_wrapper(self.k, ydomains)

        cp = daeVariable_wrapper(self.cp, domains)
        T = daeVariable_wrapper(self.T, domains)
        H = daeVariable_wrapper(self.H, domains)

        eq.Residual = k * cp * T - H


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        Di = self.Di()
        D = daeVariable_wrapper(self.D, domains)

        eq.Residual = D - Di


    def eq_Reynolds(self):

        eq = self.CreateEquation("Re", "Reynolds")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        v = daeVariable_wrapper(self.v, domains)
        Re = daeVariable_wrapper(self.Re, domains)
        D = daeVariable_wrapper(self.D, domains)
        rho = daeVariable_wrapper(self.rho, domains)
        mu = daeVariable_wrapper(self.mu, domains)

        eq.Residual = Re - D * Abs(v) * rho / mu


    def eq_fD(self):

        eq = self.CreateEquation("fDturb", "Turbulent darcy factor")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        Re = daeVariable_wrapper(self.Re, domains)
        D = daeVariable_wrapper(self.D, domains)
        fD = daeVariable_wrapper(self.fD, domains)
        ep = self.ep()

        A = (2.457 * Log(((7 / Re) ** 0.9 + 0.27 * ep / D) ** -1)) ** 16
        B = (37530 / Re) ** 16
        ff = 2 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12)

        eq.Residual = fD - 4 * ff


    def eq_lowerbound_flowrate(self):

        eq = self.CreateEquation("LowerBoundFlowrate", "Lower Bound Flowrate")

        klb = self.klb()

        if self.YDomains:
            k = self.k.array('*')
            Npipes = self.Npipes.array('*')
            eq.Residual = klb - Sum(k * Npipes)
        else:
            k = self.k()
            Npipes = self.Npipes()
            eq.Residual = klb - k * Npipes


    def eq_upperbound_flowrate(self):

        eq = self.CreateEquation("UpperBoundLowrate", "Upper Bound Flowrate")

        kub = self.kub()

        if self.YDomains:
            k = self.k.array('*')
            Npipes = self.Npipes.array('*')
            eq.Residual = kub - Sum(k * Npipes)
        else:
            k = self.k()
            Npipes = self.Npipes()
            eq.Residual = kub - k * Npipes


    def eq_lowerbound_enthalpy(self):

        eq = self.CreateEquation("LowerBoundPressure", "Lower Bound Enthalpy")

        Hlb = self.Hlb()

        if self.YDomains:
            H = self.H.array(0, '*')
            Npipes = self.Npipes.array('*')
            eq.Residual = Hlb - Sum(H * Npipes)
        else:
            H = self.H(0)
            Npipes = self.Npipes()
            eq.Residual = Hlb - H * Npipes


    def eq_upperbound_enthalpy(self):

        eq = self.CreateEquation("UpperBoundPressure", "Upper Bound Enthalpy")

        Hub = self.Hub()

        if self.YDomains:
            H = self.H.array(-1, '*')
            Npipes = self.Npipes.array('*')
            eq.Residual = Hub - Sum(H * Npipes)
        else:
            H = self.H(-1)
            Npipes = self.Npipes()
            eq.Residual = Hub - H * Npipes


    def eq_mommentum_balance(self):

        eq = self.CreateEquation("MomBal", "Momentum balance")
        ydomains = distribute_on_domains(self.YDomains, eq, eClosedClosed)
        xdomains = distribute_on_domains(self.XDomains, eq, eOpenClosed)
        domains = xdomains + ydomains

        g = self.g
        L = self.L()
        tetha = self.tetha()
        k = daeVariable_wrapper(self.k, ydomains)
        P = daeVariable_wrapper(self.P, domains)
        fD = daeVariable_wrapper(self.fD, domains)
        v = daeVariable_wrapper(self.v, domains)
        D = daeVariable_wrapper(self.D, domains)
        rho = daeVariable_wrapper(self.rho, domains)
        Re = daeVariable_wrapper(self.Re, domains)
        ep = self.ep()

        hL = 0.5 * fD * Abs(v) * v / ( D * g )

        DeltaP = g * rho * hL

        Area = 0.25 * 3.14 * D ** 2

        eq.Residual =  dt(Area * rho * v) + k * d(v, self.x, eCFDM) / L + Area * d(P, self.x, eCFDM) / L  +  Area * DeltaP +  Area * rho * g * Sin(tetha)


    def eq_heat_balance(self):

        eq = self.CreateEquation("HeatBal", "Heat balance - T")
        ydomains = distribute_on_domains(self.YDomains, eq, eClosedClosed)
        xdomains = distribute_on_domains(self.XDomains, eq, eOpenClosed)
        domains = xdomains + ydomains

        L = self.L()
        D = daeVariable_wrapper(self.D, domains)
        rho = daeVariable_wrapper(self.rho, domains)
        k = daeVariable_wrapper(self.k, ydomains)
        T = daeVariable_wrapper(self.T, domains)
        cp = daeVariable_wrapper(self.cp, domains)
        Qout = daeVariable_wrapper(self.Qout, domains)

        A = 0.25 * 3.14 * D ** 2
        eq.Residual = dt(rho * cp * A * T) + k * d( cp * T, self.x, eCFDM) / L + Qout


    def eq_water_properties(self):

        eq = self.CreateEquation("Density", "Density calculation")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        T = daeVariable_wrapper(self.T, domains)
        P = daeVariable_wrapper(self.P, domains)
        rho = daeVariable_wrapper(self.rho, domains)
        Tast = T / Constant(1 * K)
        Past = P / Constant(1 * Pa)

        rho_calc = density( Tast, Past, simplified = True) * Constant(1 * (kg ** (1))*(m ** (-3)))
        eq.Residual = rho - rho_calc


        eq = self.CreateEquation("Viscosity", "Viscosity calculation")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)
        T = daeVariable_wrapper(self.T, domains)
        P = daeVariable_wrapper(self.P, domains)
        mu = daeVariable_wrapper(self.mu, domains)
        Tast = T / Constant(1 * K)
        Past = P / Constant(1 * Pa)
        mu_calc = viscosity( Tast, Past, simplified = True)  * Constant(1 * (Pa ** (1))*(s ** (1)))
        eq.Residual = mu - mu_calc

        eq = self.CreateEquation("ThermalConductivity", "Thermal Conductivity calculation")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)
        T = daeVariable_wrapper(self.T, domains)
        P = daeVariable_wrapper(self.P, domains)
        kappa = daeVariable_wrapper(self.kappa, domains)
        Tast = T / Constant(1 * K)
        Past = P / Constant(1 * Pa)
        kappa_calc = conductivity( Tast, Past, simplified = True) * Constant(1 * (K ** (-1))*(W ** (1))*(m ** (-1)))
        eq.Residual = kappa - kappa_calc

        eq = self.CreateEquation("HeatCapacity", "Heat Capacity calculation")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)
        T = daeVariable_wrapper(self.T, domains)
        P = daeVariable_wrapper(self.P, domains)
        cp = daeVariable_wrapper(self.cp, domains)
        Tast = T / Constant(1 * K)
        Past = P / Constant(1 * Pa)
        cp_calc = heat_capacity( Tast, Past, simplified = True) * Constant(1 * (J ** (1))*(K ** (-1))*(kg ** (-1)))
        eq.Residual = cp - cp_calc


    def eq_total_he(self):
        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)
        Qout = daeVariable_wrapper(self.Qout, domains)
        eq.Residual = Qout


    def DeclareEquations(self):

        Edge.DeclareEquations(self)

        # friction factors
        self.eq_fD() # x * y

        # Nodal Specifications (P1, P2, T1)
        self.eq_pressure_boundaries() # 2 * y
        self.eq_temperature_boundaries() # y

        # Internal Diameter
        self.eq_internal_diameter() # x * y

        # v and Re
        self.eq_velocity() # x * y
        self.eq_Reynolds() # x * y

        # mu, rho, cp and kappa
        self.eq_water_properties() # 4 * x * y

        # P
        self.eq_mommentum_balance() # (x-1) * y

        # T, H, and Qout
        self.eq_entalphy() # x * y
        self.eq_heat_balance() # (x-1) * y
        self.eq_total_he()  # x * y

        # Bounds
        self.eq_upperbound_flowrate() # 1
        self.eq_lowerbound_flowrate() # 1

        self.eq_upperbound_enthalpy() # 1
        self.eq_lowerbound_enthalpy() # 1

        # total equations 4 + 12 * x * y + y (per edge)
        # total equations 2 (per node)

        # total variables 4 + 12 * x * y + y (per edge)
        # total variables 3 (per node)

        # IF 1 edge + 2 nodes                    x/y=>3/2       x/y=> 3/1

            # equations: 8 + 12 * x * y + y      82            45
            # variables: 10 + 12 * x * y + y     84            47
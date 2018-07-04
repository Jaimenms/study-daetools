__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""


from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from daetools_extended.tools import get_node_tree, execute_recursive_method, get_initialdata_from_reporter, update_initialdata, daeVariable_wrapper, distribute_on_domains

from pyUnits import m, kg, s, K, Pa, J, W, rad, day


class Biofilm(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        daeModelExtended.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)

        if not Parent:
            self.define_variables()
            self.define_parameters()

    def define_variables(self):

        # Variable types
        mass_biofilm_t = daeVariableType("mass_biofilm_t", (kg ** (1)) * (m ** (-2)), 1e-9, 1e+04, 1e-06, 1e-07)
        b_t = daeVariableType("b_t", (day ** (-1)), 1e-6, 1e+04, 0.2, 1e-06)
        rate_t = daeVariableType("rate_t", (kg ** (1)) * (m ** (-2)) * (day ** (-1)), 1e-9, 1e+04, 0.2, 1e-06)

        self.mf = daeVariable("mf", mass_biofilm_t, self, "Biofilm Mass per Area", self.Domains)
        self.Tbf = daeVariable("Tbf", temperature_t, self, "Temperature of Biofilm", self.Domains)
        self.vbf = daeVariable("vbf", velocity_t, self, "Velocity over Biofilm", self.Domains)
        self.Jp = daeVariable("Jp", rate_t, self, "Biofilm formation specific formation rate", self.Domains)
        self.b = daeVariable("b", b_t, self, "Biofilm formation specific removal rate", self.Domains)
        self.ratef = daeVariable("ratef", rate_t, self, "Biofilm formation specific removal rate", self.Domains)

    def define_parameters(self):

        self.lagt = Constant(24*3600*s)

        self.rhomf = daeParameter("rhomf",(kg ** (1)) * (m ** (-3)), self, "Density of the biofilm")
        self.mfi = daeParameter("mfi", (kg ** (1)) * (m ** (-2)), self, "Initial biofilm density")


    def eq_biofilm_Jp(self):

        eq = self.CreateEquation("BiofilmJp", "Biofilm Formation Specific Rate")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        b = daeVariable_wrapper(self.b, domains)
        Jp = daeVariable_wrapper(self.Jp, domains)
        Tbf = daeVariable_wrapper(self.Tbf, domains)
        vbf = daeVariable_wrapper(self.vbf, domains)
        vast = vbf / Constant(1 * m / s)
        Tast = Tbf / Constant(1 * K)

        Jp_27 = (-0.7006751988520601 * vast + 0.7407489411915851) / (1. + Exp(20 * ((vast - 1.0)))) * Constant(1 * kg / m ** 2 / day)

        sigmoid = 1 * Exp(0.6221 * (Tast - 315.34)) / (1 + Exp(0.6221 * (Tast - 315.34)))
        k27_mf = 0.599
        k_mf = 1.624e7 * Tast * Exp(-1000 * 13.609 / 1.987 / Tast) * (1 - sigmoid)

        Jp_model = Jp_27 * k_mf / k27_mf

        eq.Residual = Jp - Jp_model


    def eq_biofilm_rate(self):

        eq = self.CreateEquation("BiofilmRate", "Biofilm Rate")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        b = daeVariable_wrapper(self.b, domains)
        mf = daeVariable_wrapper(self.mf, domains)
        Jp = daeVariable_wrapper(self.Jp, domains)
        ratef = daeVariable_wrapper(self.ratef, domains)

        # ratef_model = (Jp - b * mf) * 0.5 * (1. + Abs(mfmax - mf) / (mfmax - mf))
        ratef_model = (Jp - b * mf)

        eq.Residual = ratef - ratef_model


    def eq_biofilm_b(self):

        eq = self.CreateEquation("Biofilmb", "Biofilm Removal Specific Rate")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        b = daeVariable_wrapper(self.b, domains)
        vbf = daeVariable_wrapper(self.vbf, domains)
        vast = vbf / Constant(1 * m / s)

        b_model = 1 / (6.708428763473813 * vast + 1.951170457848189) * Constant(1 / day)

        eq.Residual = b - b_model


    def eq_biofilm(self):

        self.IF(Time() < self.lagt, eventTolerance=1E-5)

        eq = self.CreateEquation("BiofilmOFF", "Biofilm - OFF")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        mf = daeVariable_wrapper(self.mf, domains)
        eq.Residual = mf - self.mfi()

        self.ELSE()

        eq = self.CreateEquation("BiofilmON", "Biofilm - ON")
        domains = distribute_on_domains(self.Domains, eq, eClosedClosed)

        mf = daeVariable_wrapper(self.mf, domains)
        ratef = daeVariable_wrapper(self.ratef, domains)
        eq.Residual = dt(mf) - ratef / Constant( 24*3600 * s / day)

        self.END_IF()


    def DeclareEquations(self):

        daeModelExtended.DeclareEquations(self)

        self.eq_biofilm()
        self.eq_biofilm_Jp()
        self.eq_biofilm_b()
        self.eq_biofilm_rate()
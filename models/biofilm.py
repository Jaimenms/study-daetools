__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""


from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended

from pyUnits import m, kg, s, K, Pa, J, W, rad


class Biofilm(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        daeModelExtended.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_variables(self):

        # Variable types
        mass_biofilm_t = daeVariableType("mass_biofilm_t", (kg ** (1)) * (m ** (-2)), 1e-10, 1e+04, 0.1, 1e-05)
        self.mf = daeVariable("m_f", mass_biofilm_t, self, "Mass of Biofilm", [self.x, ])


    def define_parameters(self):

        self.rhomf = daeParameter("rho_mf",(kg ** (1)) * (m ** (-3)), self, "Density of the biofilm")

        pass


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.D(x) - (self.Di() ** 2 - 4 * self.mf(x) * self.Di() / self.rhomf()) ** 0.5


    def eq_biofilm(self):

        eq = self.CreateEquation("Biofilm", "Biofilm Formation - ON")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        vast = self.v(x) / Constant(1 * m / s)
        Tast = self.T(x) / Constant(1 * K)
        Jp = (1.19e-7 - 1.14e-7 * vast) / (1.+Exp(20*((vast - 1.0)))) * Constant( 1 * kg / m **2 / s)
        b =  1 / (4.26e4  + 3.16e5 * vast) * Constant( 1 / s)
        k27_mf = 0.599
        sigmoid = 1 * Exp(0.6221 * (Tast - 315.34)) / (1 + Exp(0.6221 * (Tast - 315.34)))
        k_mf = 1.624e7 * Tast * Exp(-1000 * 13.609 / 1.987 / Tast) * (1 - sigmoid)
        #eq.Residual = dt(self.mf(x)) - ( Jp - b * self.mf(x) )
        eq.Residual = dt(self.mf(x)) - ( Jp * k27_mf / k_mf - b * self.mf(x) )


    def DeclareEquations(self):

        daeModelExtended.DeclareEquations(self)

        self.eq_biofilm()

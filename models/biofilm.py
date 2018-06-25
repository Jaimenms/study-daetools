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
        self.mf = daeVariable("mf", mass_biofilm_t, self, "Mass of Biofilm", [self.x, ])
        self.Tbf = daeVariable("Tbf", temperature_t, self, "Temperature of Biofilm", [self.x, ])


    def define_parameters(self):

        self.lagt = Constant(10 * s)

        self.rhomf = daeParameter("rhomf",(kg ** (1)) * (m ** (-3)), self, "Density of the biofilm")
        self.mfi = daeParameter("mfi", (kg ** (1)) * (m ** (-2)), self, "Initial biofilm density")


    def eq_biofilm(self):

        self.IF(Time() < self.lagt, eventTolerance=1E-5)

        eq = self.CreateEquation("BiofilmOFF", "Biofilm - OFF")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.mf(x) - self.mfi()

        self.ELSE()

        eq = self.CreateEquation("BiofilmON", "Biofilm - ON")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        vast = self.v(x) / Constant(1 * m / s)


        Tast = self.Tbf(x) / Constant(1 * K) # In case we consider temperature of constact with water
        #Tast = 0.5 * (self.Tbf(x) + self.Ti(x)) / Constant(1 * K) # In case we consider the mean temperature

        Jp = (1.19e-7 - 1.14e-7 * vast) / (1.+Exp(20*((vast - 1.0)))) * Constant( 1 * kg / m **2 / s)
        b =  1 / (4.26e4  + 3.16e5 * vast) * Constant( 1 / s)
        k27_mf = 0.599
        sigmoid = 1 * Exp(0.6221 * (Tast - 315.34)) / (1 + Exp(0.6221 * (Tast - 315.34)))
        k_mf = 1.624e7 * Tast * Exp(-1000 * 13.609 / 1.987 / Tast) * (1 - sigmoid)
        eq.Residual = dt(self.mf(x)) - ( Jp * k27_mf / k_mf - b * self.mf(x) )

        self.END_IF()



    def DeclareEquations(self):

        daeModelExtended.DeclareEquations(self)

        self.eq_biofilm()

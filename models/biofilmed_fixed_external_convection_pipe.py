__doc__ = """
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.fixed_external_convection_pipe import FixedExternalConvectionPipe
    from models.biofilm import Biofilm
except:
    from .fixed_external_convection_pipe import FixedExternalConvectionPipe
    from .biofilm import Biofilm


class BiofilmedFixedExternalConvectionPipe(Biofilm, FixedExternalConvectionPipe):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        FixedExternalConvectionPipe.__init__(self, Name, Parent=Parent, Description=Description, data=data,
                                             node_tree=node_tree)

    def define_parameters(self):
        FixedExternalConvectionPipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.D(x) - (self.Di() ** 2 - 4 * self.mf(x) * self.Di() / self.rhomf()) ** 0.5


    def eq_calculate_resistance(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Resext = 1 / (self.pi * self.Do() * self.hext())
        Resint = 1 / (self.pi * self.D(x) * self.hint(x))
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        Resfilm = Log(self.Di() / self.D()) / (2 * self.pi * self.kappa()) # Melhorar ajustando T para kappa
        eq.Residual = self.Resistance(x) - (Resint + Reswall + Resext + Resfilm)


    def eq_biofilm_temperature(self):

        eq = self.CreateEquation("Tbf", "Biofilm Temperature")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        # eq.Residual = self.Tbf(x) - 0.5 * (self.T(x) + self.Ti(x))
        #eq.Residual = self.Tbf(x) - 0.5 * (self.T(x) + self.Ti(x))
        eq.Residual = self.Qout(x) - (self.T(x) - self.Tbf(x)) * (self.pi * self.D(x) * self.hint(x))


    def define_variables(self):

        FixedExternalConvectionPipe.define_variables(self)

        Biofilm.define_variables(self)


    def define_parameters(self):

        FixedExternalConvectionPipe.define_parameters(self)

        Biofilm.define_parameters(self)


    def DeclareEquations(self):

        FixedExternalConvectionPipe.DeclareEquations(self)

        self.eq_biofilm()
        self.eq_biofilm_temperature()

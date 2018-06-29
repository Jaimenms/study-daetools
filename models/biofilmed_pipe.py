__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""


from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.pipe import Pipe
    from models.biofilm import Biofilm
except:
    from .pipe import Pipe
    from .biofilm import Biofilm


class BiofilmedPipe(Biofilm, Pipe):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        Pipe.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.D(x) - (self.Di() ** 2 - 4 * self.mf(x) * self.Di() / self.rhomf()) ** 0.5


    def eq_biofilm_temperature(self):

        eq = self.CreateEquation("Tbf", "Biofilm Temperature")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        #eq.Residual = self.Tbf(x) - 0.5 * (self.T(x) + self.Ti(x))
        eq.Residual = self.T(x) - self.Tbf(x)


    def eq_biofilm_velocity(self):
        eq = self.CreateEquation("vbf", "Biofilm Velocity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        eq.Residual = self.v(x) - self.vbf(x)


    def define_parameters(self):
        Pipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def define_variables(self):
        Pipe.define_variables(self)
        Biofilm.define_variables(self)


    def define_parameters(self) :
        Pipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def DeclareEquations(self):
        Pipe.DeclareEquations(self)

        self.eq_biofilm()
        self.eq_biofilm_Jp()
        self.eq_biofilm_b()
        self.eq_biofilm_rate()
        self.eq_biofilm_temperature()
        self.eq_biofilm_velocity()


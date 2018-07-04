__doc__ = """
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

import pandas as pd


try:
    from models.external_film_condensation_tube_arrange import ExternalFilmCondensationTubeArrange
    from models.biofilm import Biofilm
except:
    from .external_film_condensation_tube_arrange import ExternalFilmCondensationTubeArrange
    from .biofilm import Biofilm

from daetools_extended.tools import get_node_tree, execute_recursive_method, get_initialdata_from_reporter, update_initialdata

class BiofilmedExternalFilmCondensationTubeArrange(Biofilm, ExternalFilmCondensationTubeArrange):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        ExternalFilmCondensationTubeArrange.__init__(self, Name, Parent=Parent, Description=Description, data=data,
                                             node_tree=node_tree)

    def define_parameters(self):
        ExternalFilmCondensationTubeArrange.define_parameters(self)
        Biofilm.define_parameters(self)


    def eq_internal_diameter(self):

        eq = self.CreateEquation("D", "D_internal_flow_diameter")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        y = eq.DistributeOnDomain(self.y, eClosedClosed)
        eq.Residual = self.D(x,y) - (self.Di() ** 2 - 4 * self.mf(x,y) * self.Di() / self.rhomf()) ** 0.5


    def eq_calculate_resistance(self):

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        y = eq.DistributeOnDomain(self.y, eClosedClosed)
        Resext = 1 / (self.pi * self.Do() * self.hext(x,y))
        Resint = 1 / (self.pi * self.D(x,y) * self.hint(x,y))
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        Resfilm = Log(self.Di() / self.D(x,y)) / (2 * self.pi * self.kappa(x,y)) # Melhorar ajustando T para kappa
        eq.Residual = self.Resistance(x,y) - (Resint + Reswall + Resext + Resfilm)


    def eq_biofilm_temperature(self):

        eq = self.CreateEquation("Tbf", "Biofilm Temperature")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        y = eq.DistributeOnDomain(self.y, eClosedClosed)
        #eq.Residual = self.Qout(x) - (self.T(x) - self.Tbf(x)) * (self.pi * self.D(x) * self.hint(x))
        eq.Residual = self.T(x, y) - self.Tbf(x, y)

    def eq_biofilm_velocity(self):

        eq = self.CreateEquation("vbf", "Biofilm Velocity")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        y = eq.DistributeOnDomain(self.y, eClosedClosed)
        eq.Residual = self.v(x, y) - self.vbf(x, y)


    def define_variables(self):
        ExternalFilmCondensationTubeArrange.define_variables(self)

        Biofilm.define_variables(self)


    def define_parameters(self):
        ExternalFilmCondensationTubeArrange.define_parameters(self)

        Biofilm.define_parameters(self)


    def DeclareEquations(self):

        ExternalFilmCondensationTubeArrange.DeclareEquations(self)
        Biofilm.DeclareEquations(self)

        self.eq_biofilm_temperature()
        self.eq_biofilm_velocity()
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

from water_properties import density, viscosity, conductivity, heat_capacity


class BiofilmedFixedExternalConvectionPipe(Biofilm, FixedExternalConvectionPipe):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        FixedExternalConvectionPipe.__init__(self, Name, Parent=Parent, Description=Description, data=data,
                                             node_tree=node_tree)

    def define_parameters(self):
        FixedExternalConvectionPipe.define_parameters(self)
        Biofilm.define_parameters(self)

    def eq_total_he(self):
        print("PASSOU AQUI")

        eq = self.CreateEquation("TotalHeat", "Heat balance - Qout")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)

        Resext = 1 / (2 * self.pi * self.Di() * self.hext())
        Resint = 1 / (2 * self.pi * self.D(x) * self.hint(x))
        Reswall = Log(self.Do() / self.Di()) / (2 * self.pi * self.kwall())
        Tm = 0.5 * (self.T(x) + self.Tf(x))
        kappaf = conductivity(Tm / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified=True)
        Resf = Log(self.Di() / self.D()) / (2 * self.pi * kappaf)

        eq.Residual = self.Qout(x) - (self.T(x) - self.Text()) / (Resint + Reswall + Resext + Resf)


    def eq_calculate_Tf(self):

        eq = self.CreateEquation("FilmWall", "Heat balance - biofilm wall")
        x = eq.DistributeOnDomain(self.x, eClosedClosed)
        Tm = 0.5 * (self.T(x) + self.Tf(x))
        kappaf = conductivity(Tm / Constant(1 * K), self.P(x) / Constant(1 * Pa), simplified=True)
        eq.Residual = self.Qout(x) * (self.Di() - self.D(x)) - (self.Tf(x) - self.Ti(x)) * self.pi * self.D(x) * kappaf * Constant(
            1 * (K ** (-1)) * (W ** (1)) * (m ** (-1)))

    def define_variables(self):
        FixedExternalConvectionPipe.define_variables(self)
        Biofilm.define_variables(self)

        # State variables
        water_temperature_t = daeVariableType("temperature_t", (K ** (1)), 273.0, 400.0, 300.0, 0.01)
        self.Tf = daeVariable("Tf", water_temperature_t, self, "Film Temperature", [self.x, ])

    def define_parameters(self):
        FixedExternalConvectionPipe.define_parameters(self)
        Biofilm.define_parameters(self)

    def DeclareEquations(self):
        FixedExternalConvectionPipe.DeclareEquations(self)

        Biofilm.eq_biofilm(self)

        self.eq_calculate_Tf()

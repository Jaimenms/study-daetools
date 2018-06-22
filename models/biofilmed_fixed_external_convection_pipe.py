__doc__="""
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
        FixedExternalConvectionPipe.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_parameters(self):
        FixedExternalConvectionPipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def define_variables(self):
        FixedExternalConvectionPipe.define_variables(self)
        Biofilm.define_variables(self)


    def define_parameters(self) :
        FixedExternalConvectionPipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def DeclareEquations(self):

        FixedExternalConvectionPipe.DeclareEquations(self)

        Biofilm.eq_biofilm(self)

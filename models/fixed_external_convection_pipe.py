__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""


from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.isolated_pipe import IsolatedPipe
    from models.fixed_external_convection import FixedExternalConvection
except:
    from .isolated_pipe import IsolatedPipe
    from .fixed_external_convection import FixedExternalConvection


class FixedExternalConvectionPipe(FixedExternalConvection, IsolatedPipe):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        IsolatedPipe.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_parameters(self):

        IsolatedPipe.define_parameters(self)
        FixedExternalConvection.define_parameters(self)


    def define_variables(self):

        IsolatedPipe.define_variables(self)
        FixedExternalConvection.define_variables(self)


    def DeclareEquations(self):

        IsolatedPipe.DeclareEquations(self)
        FixedExternalConvection.eq_total_he(self)
        FixedExternalConvection.eq_wall_he(self)



__doc__="""
DaetTools model that describes the behavior of a water flowing in a pipe with the effect of biofim formation.
"""


from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.isolated_pipe import IsolatedPipe
    from models.biofilm import Biofilm
except:
    from .isolated_pipe import IsolatedPipe
    from .biofilm import Biofilm


class BiofilmedIsolatedPipe(Biofilm, IsolatedPipe):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):

        IsolatedPipe.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def define_parameters(self):
        IsolatedPipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def define_variables(self):
        IsolatedPipe.define_variables(self)
        Biofilm.define_variables(self)


    def define_parameters(self) :
        IsolatedPipe.define_parameters(self)
        Biofilm.define_parameters(self)


    def DeclareEquations(self):
        IsolatedPipe.DeclareEquations(self)
        Biofilm.eq_biofilm(self)


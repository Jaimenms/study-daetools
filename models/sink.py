__doc__="""
DaetTools model that describes the behavior of a water sink node without dynamics effects. It has tree state variables:

* P : pressure in Pa
* T : temperature in K
* w : inlet nodal mass flowrate in kg/s

"""

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

try:
    from models.node import Node
except:
    from .node import Node

from water_properties import heat_capacity


class Sink(Node):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        """
        Model for the sink
        :param Name: name of the model
        :param Parent: parent model
        :param Description: description of the model
        :param data: parameters and other required data
        """

        # Instantiate the Node
        Node.__init__(self, Name, Parent=Parent, Description=Description, data=data, node_tree=node_tree)


    def eq_mass_balance(self):
        """
        This method writes the mass balance to the correspondent node instance
        :return:
        """

        # Starting with the external mass flow rate
        residual_aux = -self.w()

        # Mass to the node inlet
        for edge_name in self.get_inlet():
            residual_aux += self.Parent.submodels[edge_name].k()
            print("++ edge {0} upstream".format(edge_name))

        # Mass to the node outlet
        for edge_name in self.get_outlet():
            residual_aux -= self.Parent.submodels[edge_name].k()
            print("++ edge {0} downstream".format(edge_name))

        # Instantiate equation NMB
        eq = self.CreateEquation("NMB_nodal_mass_balance")
        eq.Residual = residual_aux

        print("+ sink mass balance")


    def eq_energy_balance(self):
        """
        This method writes the mass balance to the correspondent node instance
        :return:
        """

        # Starting with the external mass flow rate
        cp_nodal = heat_capacity( self.T() / Constant(1 * K), self.P() / Constant(1 * Pa), simplified = True)
        residual_aux = -self.w() * self.T() * cp_nodal * Constant(1 * (J ** (1))*(K ** (-1))*(kg ** (-1)))

        # Mass to the node inlet
        for edge_name in self.get_inlet():
            ind_edge_out = self.Parent.submodels[edge_name].domain_len - 1
            residual_aux += self.Parent.submodels[edge_name].H(ind_edge_out)

        # Mass to the node outlet
        for edge_name in self.get_outlet():
            # TODO - Consider cp in node and not external
            residual_aux -= self.Parent.edges[edge_name].H(0)

        eq = self.CreateEquation("NEB_nodal_energy_balance_2")
        eq.Residual = residual_aux
        print("+ sink energy_balance 2")


    def DeclareEquations(self):
        """
        This Method is called by the DaeTools. Here is where all the equations are defined
        :return:
        """

        Node.DeclareEquations(self)

        print("Reading Sink Equations")
        self.eq_mass_balance()
        self.eq_energy_balance()

__doc__="""
This is an extension of daeModel class that permits the  automatic construction of network models based on a
data dictionary. The models can assume two different kinds:

* node: if it describes a node in a digraph representation
* edge: if it describes a edge in a digraph representation
"""

# TODO : Improve documentation

from daetools.pyDAE import *
from pyUnits import m, kg, s, K, Pa, J, W, rad

import numpy as np

from .tools import get_module_class_from_data


class daeModelExtended(daeModel):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        """
        Model for the network
        :param Name: name of the model
        :param Parent: parent model
        :param Description: description of the model
        :param data: parameters and other required data
        """

        # If data has name attribute, it overwrites the Name variable
        if "name" in data:
            Name = data['name']

        # If data has description attribute, it overwrites the Name variable
        if "description" in data:
            Description = data['description']

        # Initialize the daeModel class
        daeModel.__init__(self, Name, Parent, Description)


        # If it is nodal
        if 'kind' in data and self.check_if_model_is_nodal(data['kind']):
            # Set list of edges that sends fluid to a node with name equal to submodel_name
            data = self.set_inlet_outlet(data['kind'], Name, data, node_tree, position='inlet')

            # Set list of edges that takes fluid from a node with name equal to submodel_name
            data = self.set_inlet_outlet(data['kind'], Name, data, node_tree, position='outlet')

        # Collect the data dictionary
        self.data = data

        # Collect Parent
        self.Parent = Parent

        # Read the submodels tree
        self.instantiate_submodels(node_tree)


    def instantiate_submodels(self, node_tree):
        """
        This method reads the data dict to instantiate the submodels if data contains the element submodels as a dictionary
        :return:
        """

        self.submodels = {}

        # If there is no submodels, return
        if "submodels" not in self.data:
            return

        # Test if submodels structure is readable
        # TODO - Include here a better error handling
        if not type(self.data['submodels']) == dict:
            print("It was not possible to collect submodels information")
            return

        # Instantiate each submodel
        for submodel_name, submodel_data in self.data["submodels"].items():
            print("Instantiating {0}".format(submodel_name))
            self.instantiate_submodel(submodel_name, submodel_data, node_tree)



    def instantiate_submodel(self, submodel_name, submodel_data, node_tree):
        """
        This method instantiate each submodel, based on the model dictionary
        :return:
        """

        class_ = get_module_class_from_data(submodel_data)

        # get relevant data
        module_name = submodel_data['module']
        class_name = submodel_data['class']

        # Instantiate the submodel class
        print("with class {0}".format(class_name))
        self.submodels[submodel_name] = class_(
            submodel_name,
            Parent=self,
            Description="Submodel {0} with {1}.{2}".format(submodel_name, module_name, class_name),
            data=submodel_data,
            node_tree=node_tree,
        )


    def set_inlet_outlet(self, kind, name, submodel_data, node_tree, position='inlet'):
        """
        Sets the list of edges that sends/gets fluid to/from the node
        :param kind:
        :param name:
        :return:
        """

        if self.check_if_submodel_in_tree(kind, name, node_tree):

            print("Setting position {0} to the node {1}".format(position,name), node_tree[name][position])

            submodel_data[position] = node_tree[name][position]

        return submodel_data


    def check_if_submodel_in_tree(self, submodel_kind, submodel_name, node_tree):
        """
        Check if submodel is a node that is contained in the node_tree
        :param submodel_name:
        :param submodel_kind:
        :return:
        """

        # Check if this node is in the tree
        if submodel_name in node_tree:
            return True
        else:
            return False


    def check_if_model_is_nodal(self,submodel_kind):
        """

        :param submodel_kind: check if the kind is nodal
        :return:
        """

        NODAL_KINDS = ('node',)

        if submodel_kind in NODAL_KINDS:
            return True
        else:
            return False


    def setup_domains(self):
        """
        Setup domains according to data dictionary structure
        :return:
        """

        if 'domains' in self.data:
            self.x.CreateStructuredGrid(self.domain_len-1, 0.0, 1.0)


    def setup_active_states(self):
        """
        Setup active states according to data dictionary structure
        :return:
        """

        if 'states' in self.data:
            for name, value in self.data['states'].items():
                getattr(self, name).ActiveState = value


    def setup_variables(self):
        """
        Setup variables according to data dictionary structure
        :return:
        """

        # Setting the specifications
        if 'specifications' in self.data:
            for name, value in self.data['specifications'].items():
                n = getattr(self, name).NumberOfPoints
                print(name,n)
                if n == 1:
                    getattr(self, name).AssignValue(value)
                else:
                    for i in range(n):
                        getattr(self, name).AssignValue(i, value)

        if 'initial_conditions' in self.data:
            for name, value in self.data['initial_conditions'].items():
                n = getattr(self, name).NumberOfPoints
                if n == 1:
                    getattr(self, name).SetInitialCondition(value)
                else:
                    for i in range(n):
                        getattr(self, name).SetInitialCondition(i, value)


    def setup_parameters(self):
        """
        Setup parameters according to data dictionary structure
        :return:
        """
        # Setting the parameters
        if 'parameters' in self.data:
            for name, value in self.data['parameters'].items():
                getattr(self, name).SetValue(value)


    def setup_initial_guess(self):
        """
        Setup initial guesses according to data dictionary structure
        :return:
        """

        # Setting the parameters
        if 'initial_guess' in self.data:
            for name, value in self.data['initial_guess'].items():
                n = getattr(self, name).NumberOfPoints
                if type(value) == list and len(value) == 2:
                    i = 0
                    for valuei in np.linspace(value[0],value[1],n):
                        getattr(self, name).SetInitialGuess(i, valuei)
                        i += 1

                elif type(value) == list and len(value) > 2:
                    i = 0
                    for valuei in value:
                        getattr(self, name).SetInitialGuess(i, valuei)
                        i += 1
                else:
                    if n == 1:
                        getattr(self, name).SetInitialGuess(value)
                    else:
                        for i in range(0,n):
                            getattr(self, name).SetInitialGuess(i, value)


    def get_inlet(self):
        """
        Get the list of all edges that gives fluid to the node
        :return:
        """
        return self.data['inlet']


    def get_outlet(self):
        """
        Get the list of all edges that takes fluid from the node
        :return:
        """

        return self.data['outlet']


    def get_from(self):
        """
        Get the node from which the fluid enters the edge
        :return:
        """

        return self.data['from']


    def get_to(self):
        """
        Get the node which the fluid goes
        :return:
        """

        return self.data['to']


    def get_node(self, position):
        """
        Get the node acoording to the position of the edge (from or to)
        :param position:
        :return:
        """

        if position == 'from':
            return self.get_from()
        elif position == 'to':
            return self.get_to()
        else:
            return None


    def DeclareEquations(self):
        """
        DeclareEquations mandatory for DaeTools
        :return:
        """

        daeModel.DeclareEquations(self)

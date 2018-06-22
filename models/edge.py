__doc__="""
DaetTools base edge model to be extended by pipe and others:
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended
from pyUnits import m, kg, s, K, Pa, mol, J, W, rad


class Edge(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}, node_tree={}):
        """
        Model for the edge elements
        :param Name: name of the model
        :param Parent: parent model
        :param Description: description of the model
        :param data: parameters and other required data
        """

        daeModelExtended.__init__(self, Name, Parent, Description, data=data, node_tree=node_tree)

        # Domains
        self.define_domains()

        # Defining Constants
        self.define_constants()

        # Defining Parameters
        self.define_parameters()

        # Defining Variables
        self.define_variables()


    def define_domains(self):

        pass


    def define_constants(self):

        pass


    def define_parameters(self):

        pass


    def define_variables(self):

        pass


    def eq_nodal_boundaries(self, name='P', description='pressure', position='from'):

        nodename = self.get_node(position)
        if nodename:
            if position == 'from':
                i = 0
            else:
                i = getattr(self, name).NumberOfPoints - 1
            # Instantiate equation PCB1
            eq = self.CreateEquation("BC_{0}_{1}_{2}_{3}".format(name, description, position, nodename))
            eq.Residual = getattr(self, name)(i) - getattr(self.Parent.submodels[nodename],name)()
            print("+ edge BC_{0}_{1}_{2}_{3}".format(name, description, position, nodename))


    def eq_generic_boundaries(self, name='P', description='pressure', position='bc-ini', unit=(Pa)):

        if position in self.data:
            for namei, value in self.data[position].items():
                if namei == name:
                    if position == 'bc-ini':
                        i = 0
                    else:
                        i = getattr(self, name).NumberOfPoints - 1
                    eq = self.CreateEquation("BC_{0}_{1}_{2}".format(name,description,i))
                    eq.Residual = getattr(self, name)(i) - Constant(value * unit)
                    print("+ edge " + position)


    def eq_pressure_boundaries(self):
        """
        This method writes the pressure balance to the correspondent node instance
        :return:
        """

        self.eq_generic_boundaries(name='P', description='pressure', position='bc-ini', unit=(Pa))
        self.eq_generic_boundaries(name='P', description='pressure', position='bc-fini', unit=(Pa))
        self.eq_nodal_boundaries(name='P', description='pressure', position='from')
        self.eq_nodal_boundaries(name='P', description='pressure', position='to')


    def eq_temperature_boundaries(self):
        """
        This method writes the pressure balance to the correspondent node instance
        :return:
        """

        self.eq_generic_boundaries(name='T', description='temperature', position='bc-ini', unit=(K))
        self.eq_generic_boundaries(name='T', description='temperature', position='bc-fini', unit=(K))
        self.eq_nodal_boundaries(name='T', description='temperature', position='from')


    def DeclareEquations(self):
        """
        This Methos is called by the DaeTools. Here is where all the equations are defined
        :return:
        """

        daeModelExtended.DeclareEquations(self)

        print("Reading Edge Equations")

from daetools.pyDAE import *

"""
- system.py
  - node.py
  - edge.py
    - cooling_tower.py
    - pump.py
    - pipe.py
      - fluid.py
	  - flow.py
	  - biofilm.py
	  - wall.py
	  - environment.py
	- condenser.py
	  - flow.py
	  - biofilm.py
	  - wall.py
	  - shell.py 


heat transfer:

o---)---)---) env
T0  T1  T2  T3  Tenv

res = q - f0(T0,T1)
res = q - f1(T1,T2)
res = q - f2(T2,T3)
res = q - f3(T3,Tenv)


Variables: T1, T2, T3, q
Parameter: Tenv
Fixed variable: T0

"""

from pyUnits import m, kg, s, K, Pa, mol, J, W

from models.pipe import Pipe
from models.condenser import Condenser

network_data = dict()
network_data.nodes = {
    'node-A': {
        'id': 0,
        'P': 1.0,
    },
    'node-B': {
        'id': 1,
        'w', 0.0,
    },
    'node-C': {
        'id': 2,
        'w', 0.0
    },
    'node-D': {
        'id': 3,
        'P', 0.995
    }
}

network_data.edges = {
    'pipe-01': {
        'from': 'node-A',
        'to': 'node-B',
        'module': 'pipe',
        'class': 'Pipe',
    },
    'condenser': {
        'from': 'node-B',
        'to': 'node-C',
        'module': 'condenser',
        'class': 'Condenser',
    },
    'pipe-02': {
        'from': 'node-C',
        'to': 'node-D',
        'module': 'pipe',
        'class': 'Pipe',
    }
}


def calculate_network_interconnection(network_data):

    for node_name, node_data in network_data.nodes.items():

        if 'inlet' not in node_data:
            node_data['inlet'] = []

        if 'outlet' not in node_data:
            node_data['outlet'] = []

        for edge_name, edge_data in network_data.edge.items():

            if edge_data['from'] == node_name:
                node_data['outlet'].append(edge_name)

            if edge_data['to'] == node_name:
                node_data['inlet'].append(edge_name)


class Edges(daeModel):


    def __init__(self, Name, Parent=None, Description="", data={}):

        import importlib

        daeModel.__init__(self, Name, Parent, Description)

        # Defining Variables
        self.Pin = daeVariable("P_in", no_t, self)
        self.Pin.Description = "Pressure in edge inlet"

        self.Pout = daeVariable("P_out", no_t, self)
        self.Pout.Description = "Pressure in edge outlet"

        # TODO - Verificar se isso está funcionando
        module_name = data['module']
        class_name = data['class']
        module = importlib.import_module(module_name)
        class_ = getattr(module, class_name)
        self[class_name] = class_()


class Nodes(daeModel):

    def __init__(self, Name, Parent=None, Description="", data={}):

        daeModel.__init__(self, Name, Parent, Description)

        self.N = daeDomain("N", self, unit(), "Number of nodes")

        # Defining Variables
        self.P = daeVariable("P", no_t, self)
        self.P.Description = "Nodal Pressure"
        self.P.DistributeOnDomain(self.N)

        # Defining Variables
        self.w = daeVariable("w", mass_flowrate_t, self)
        self.w.Description = "External mass flowrate"
        self.w.DistributeOnDomain(self.N)



class Network(daeModel):

    def __init__(self, Name, Parent=None, Description="", data={}):

        daeModel.__init__(self, Name, Parent, Description)


        self.data = calculate_network_interconnection(self.data)


        # Reading Edge Models
        self.edges = {}
        for edge_name, edge_data in self.data.edges.items():
            self.edges[edge_name] = Edges(
                edge_name,
                Description="Edge {0} model".format(edge_name),
                data=edge_data,
            )

        # Reading Node Models
        self.nodes = Nodes('nodes', Description="Nodes model", data=self.data.nodes)

    def DeclareEquations(self):

        daeModel.DeclareEquations(self)

        i = 0
        for edge_name, edge_data in self.data.edges.items():

            node_name_from = edge_data['from']
            node_name_to = edge_data['to']

            node_id_from = node_name_from['id']
            node_id_to = node_name_to['id']


            eq = self.CreateEquation("Edge_Pin_{0}_to_Node".format(i))
            eq.Residual = self.edges[edge_name].Pin() - self.node.P(node_id_from)

            eq = self.CreateEquation("Edge_Pout_{0}_to_Node".format(i))
            eq.Residual = self.edges[edge_name].Pout() - self.node.P(node_id_to)


        for node_name, node_data in self.data.nodes.items():

            node_id = node_data['id']

            inlet_edges = node_data['inlet']
            outlet_edges = node_data['outlet']

            eq = self.CreateEquation("Node_{0}_Mass_Balance".format(i))
            eq.Residual = self.node.w(node_id)
            for inlet_edge in inlet_edges:
                eq.Residual += self.edges[inlet_edge].kout()
            for outlet_edge in outlet_edges:
                eq.Residual += -self.edges[outlet_edge].kin()


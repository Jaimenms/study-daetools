__doc__ = """
This is the network model. A network is digraph composed by two entities: nodes and edges. An edge is an interconnection
between two nodes. An node is an interconnection between various edges and can be used as an inlet or an outlet of mass.
A node can be an equipment like a pressure vessel, a tank, a CSTR reator or any other equipment that can be modeled with
concentrated various. An edge can be a pipeline, a heat exchanger or a plug flow reator that can be modeled with distri-
buted variables.

The equations of the network model consist of the mass balance at nodes and pressure equalization at nodes. (obs.: this
equations could be included in the node model but they are easier coded here)
"""

from daetools.pyDAE import *
from daetools_extended.daemodel_extended import daeModelExtended
from daetools_extended.tools import get_node_tree, execute_recursive_method

import pandas as pd

class Network(daeModelExtended):

    def __init__(self, Name, Parent=None, Description="", data={}):
        """
        Model for the network
        :param Name: name of the model
        :param Parent: parent model
        :param Description: description of the model
        :param data: parameters and other required data
        """

        # Calculate the node tree
        node_tree = get_node_tree(Name, data)

        # Instantiate the Extended daeModel Class
        daeModelExtended.__init__(self, Name, Parent, Description, data=data, node_tree=node_tree)


    def DeclareEquations(self):
        """
        This Methos is called by the DaeTools. Here is where all the equations are defined
        :return:
        """

        # Reading equations
        daeModelExtended.DeclareEquations(self)


class NetworkSimulation(daeSimulation):

    def __init__(self, Name, Parent=None, Description="", data={}, set_reporting = True, reporting_interval = 100, time_horizon = 0):

        daeSimulation.__init__(self)

        self.m = Network(Name, Parent=Parent, Description=Description, data=data)

        self.m.SetReportingOn(set_reporting)

        if set_reporting:
            self.ReportingInterval = 100

        if time_horizon > 0:
            self.TimeHorizon = time_horizon


    def SetUpVariables(self):

        self.InitialConditionMode = eQuasiSteadyState

        execute_recursive_method(self.m,'setup_active_states')
        execute_recursive_method(self.m,'setup_variables')
        execute_recursive_method(self.m,'setup_initial_guess')


    def SetUpParametersAndDomains(self):

        execute_recursive_method(self.m,'setup_domains')
        execute_recursive_method(self.m,'setup_parameters')


def simulate(data={}):

    from daetools.pyDAE.data_reporters import daePandasDataReporter

    cfg = daeGetConfig()
    #print(cfg)
    cfg.SetBoolean('daetools.activity.printHeader', False)
    #cfg.SetBoolean('daetools.core.printInfo', True)
    #cfg.SetBoolean('daetools.IDAS.printInfo', True)
    #cfg.SetInteger('daetools.IDAS.MaxNumStepsIC', 500)
    #cfg.SetInteger('daetools.IDAS.MaxNumJacsIC', 500)
    #cfg.SetInteger('daetools.IDAS.MaxNumItersIC', 500)
    #cfg.SetBoolean('daetools.IDAS.LineSearchOffIC', True)
    cfg.SetString('daetools.core.equations.evaluationMode', 'computeStack_OpenMP')
    #cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')

    # Instantiate
    simulation = NetworkSimulation("", data=data, set_reporting = True, reporting_interval = 100, time_horizon = 10)

    # Configurate
    datareporter = daeDelegateDataReporter()

    dr1 = daeDataReporterLocal()
    dr2 = daePandasDataReporter()

    datareporter.AddDataReporter(dr1)
    datareporter.AddDataReporter(dr2)

    solver = daeIDAS()
    solver.RelativeTolerance = 1e-5
    log = daePythonStdOutLog()

    # Save the model report and the runtime model report
    # simulation.m.SaveModelReport("test.xml")
    # simulation.m.SaveRuntimeModelReport("test-rt.xml")

    # Initialize
    simulation.Initialize(solver, datareporter, log)
    print("++ Number of variables: {0}".format(simulation.TotalNumberOfVariables))
    print("++ Number of equations: {0}".format(simulation.NumberOfEquations))

    #out=daeCodeGeneratorAnalyzer()
    #out.analyzeSimulation(simulation)
    #pprint(out.runtimeInformation)

    # Solve at time = 0
    simulation.SolveInitial()

    # Run
    simulation.Run()

    # Clean up
    simulation.Finalize()

    return (simulation, dr1, dr2)


if __name__ == "__main__":

    import examples.network_examples as ex

    data = ex.case_01()

    simulation1, dr1_1, dr2_1 = simulate(data=data)

    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
       print(dr2_1.data_frame)
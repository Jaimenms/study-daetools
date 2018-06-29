import pytest
import pandas as pd

from daetools.pyDAE import *

from daetools_extended.tools import update_initialdata, get_initialdata_from_reporter
from daetools_extended.tools import get_node_tree, get_initialdata_from_reporter, update_initialdata
from daetools_extended.daesimulation_extended import daeSimulationExtended

from daetools.pyDAE.data_reporters import daePandasDataReporter


def get_testdata(case = "all"):
    """
    Get all test data from the network_examples module
    :return:
    """

    import examples.biofilm_examples as amodule

    testdata=list()
    for function_str in dir(amodule):

        if function_str[0:5] != "case_":
            continue

        if type(case) == list and function_str in case:
            function_ = getattr(amodule, function_str)
            testdata.append(function_())

        elif case == "all" or function_str == case:
            function_ = getattr(amodule, function_str)
            testdata.append(function_())

    return testdata



def simulate(data={}):

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

    # Calculate the node tree
    Name = data['name']
    node_tree = get_node_tree(Name, data)
    print(data)

    # Instantiate
    simulation = daeSimulationExtended(Name, data=data, node_tree=node_tree, set_reporting = True, reporting_interval = 3600*24, time_horizon = 3600*24*14)

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


@pytest.mark.parametrize("data", get_testdata(case="all"))
def test_biofilm(data):
    """
    Check if the reading node function can collect the correct data
    :return:
    """

    if not data:
        print("Data is empty")
        assert False

    simulation1, dr1_1, dr2_1 = simulate(data=data)

    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
       print(dr2_1.data_frame)

    assert False


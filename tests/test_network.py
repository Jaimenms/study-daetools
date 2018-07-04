import pytest
import pandas as pd

from daetools.pyDAE import *
from daetools_extended.daesimulation_extended import daeSimulationExtended
from daetools_extended.tools import get_node_tree, get_initialdata_from_reporter, update_initialdata

from daetools.pyDAE.data_reporters import daePandasDataReporter


def get_testdata(case = "all"):
    """
    Get all test data from the network_examples module
    :return:
    """

    import examples.network_examples as amodule

    testdata=list()
    for function_str in dir(amodule):

        function_ = None

        if function_str[0:5] != "case_":
            continue

        if type(case) == list and function_str in case:
            function_ = getattr(amodule, function_str)

        elif case == "all" or function_str == case:
            function_ = getattr(amodule, function_str)

        str = function_()
        testdata.append(str)
        print(function_str, str)

    return testdata


def simulate(data={}):

    cfg = daeGetConfig()
    cfg.SetBoolean('daetools.activity.printHeader', False)
    cfg.SetFloat('daetools.IDAS.MaxStep', 100.0)
    cfg.SetFloat('daetools.IDAS.relativeTolerance', 1e-3)
    cfg.GetInteger('daetools.IDAS.MaxNumSteps', 100000)

    # Calculate the node tree
    Name = data['name']
    node_tree = get_node_tree(Name, data)

    # Instantiate
    simulation = daeSimulationExtended(Name, data=data, node_tree=node_tree, set_reporting = True, reporting_interval = 3600*24, time_horizon = 3600*24)

    # Configurate
    datareporter = daeDelegateDataReporter()

    dr1 = daeDataReporterLocal()
    dr2 = daePandasDataReporter()

    datareporter.AddDataReporter(dr1)
    datareporter.AddDataReporter(dr2)

    solver = daeIDAS()
    #solver.RelativeTolerance = 1e-5
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
#@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_external_film_cond_pipe"))
def test_simulation(data):
    """
    Check if the reading node function can collect the correct data
    :return:
    """

    if not data:
        print("Data is empty")
        assert False

    simulation1, dr1_1, dr2_1 = simulate(data=data)

    previous_output = get_initialdata_from_reporter(dr1_1)

    new_data = update_initialdata(data['name'], previous_output, data)

    simulation2, dr1_2, dr2_2 = simulate(data=new_data)

    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
       print(dr2_2.data_frame)

    assert dr2_2.data_frame.loc['pipe_01.P','Values'][0][0] < 450000.


def test_simulation_advanced():
    """
    Check if the reading node function can collect the correct data
    :return:
    """
    import examples.network_examples as ex

    data1 = ex.case_external_film_condensation_pipe()

    data2 = ex.case_external_film_condensation_pipe()

    print("*******Simulation 1 started")
    simulation1, dr1_1, dr2_1 = simulate(data=data1)
    print("*******Simulation 1 concluded")

#    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
#       print(dr2_1.data_frame)

    previous_output = get_initialdata_from_reporter(dr1_1)

    new_data2 = update_initialdata(data1['name'], previous_output, data2)

    print("*******Simulation 2 started")
    simulation2, dr1_2, dr2_2 = simulate(data=new_data2)
    print("*******Simulation 2 concluded")

    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
       print(dr2_2.data_frame)

    assert dr2_2.data_frame.loc['pipe_01.P','Values'][0][0] < 450000.
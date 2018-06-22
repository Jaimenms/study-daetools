import pytest
import pandas as pd

from daetools.pyDAE import *

from models.network import *

from daetools_extended.tools import update_initialdata, get_initialdata_from_reporter


def get_testdata(case = "all"):
    """
    Get all test data from the network_examples module
    :return:
    """

    import examples.network_examples as amodule

    testdata=list()
    for function_str in dir(amodule):
        if function_str[0:5] != "case_":
            continue

        if case == "all" or function_str == case:
            function_ = getattr(amodule, function_str)
            testdata.append(function_())

    return testdata


def test_network_creation():
    """
    This first test is just to check if the basic infrastructure is up
    """

    net = Network("Creation_test")

    assert True


#@pytest.mark.skip(reason="no way of currently testing this")
#@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_fixed_external_convection_pipe"))
@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_isolated_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_isolated_pipe"))
def test_simulation(data):
    """
    Check if the reading node function can collect the correct data
    :return:
    """

    simulation1, dr1_1, dr2_1 = simulate(data=data)

    previous_output = get_initialdata_from_reporter(dr1_1)

    new_data = update_initialdata(data['name'], previous_output, data)

    new_data['submodels']['pipe_01']['states']['stnWaterPropertiesSetup'] = 'Variable'

    simulation2, dr1_2, dr2_2 = simulate(data=new_data)

    with pd.option_context('display.max_rows', None, 'display.max_columns', 20):
       print(dr2_2.data_frame)

    assert dr2_2.data_frame.loc['pipe_01.P','Values'][0][0] > 990000.
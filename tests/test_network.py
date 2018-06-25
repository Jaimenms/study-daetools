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

        if type(case) == list and function_str in case:
            function_ = getattr(amodule, function_str)
            testdata.append(function_())

        elif case == "all" or function_str == case:
            function_ = getattr(amodule, function_str)
            testdata.append(function_())

    return testdata


#@pytest.mark.skip(reason="no way of currently testing this")
#@pytest.mark.parametrize("data", get_testdata(case=["case_pipe","case_biofilmed_pipe","case_fixed_external_convection_pipe","case_biofilmed_fixed_external_convection_pipe"]))
#@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_external_film_cond_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_biofilmed_fixed_external_convection_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_external_film_condensation_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_fixed_external_convection_pipe"))
#@pytest.mark.parametrize("data", get_testdata(case="case_pipe"))
@pytest.mark.parametrize("data", get_testdata(case="all"))
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

    print(dr2_2.data_frame.loc['pipe_01.P','Values'][0])

    assert dr2_2.data_frame.loc['pipe_01.P','Values'][0][0] > 990000.
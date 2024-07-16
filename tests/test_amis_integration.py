import pytest
from sch_simulation.amis_integration.amis_integration import extract_relevant_results, returnYearlyPrevalenceEstimate, FixedParameters, run_model_with_parameters
import pandas as pd
from pandas import testing as pdt
from numpy import testing as npt

from sch_simulation.helsim_FUNC_KK.file_parsing import parse_coverage_input

example_parameters = FixedParameters(
        # the higher the value of N, the more consistent the results will be
        # though the longer the simulation will take
        number_hosts = 10,
        # no intervention
        coverage_file_name = "mansoni_coverage_scenario_0.xlsx",
        demography_name = "UgandaRural",
        # cset the survey type to Kato Katz with duplicate slide
        survey_type = "KK2",
        parameter_file_name = "mansoni_params.txt",
        coverage_text_file_storage_name = "Man_MDA_vacc.txt",
        # the following number dictates the number of events (e.g. worm deaths)
        # we allow to happen before updating other parts of the model
        # the higher this number the faster the simulation
        # (though there is a check so that there can't be too many events at once)
        # the higher the number the greater the potential for
        # errors in the model accruing.
        # 5 is a reasonable level of compromise for speed and errors, but using
        # a lower value such as 3 is also quite good
        min_multiplier = 5
    )

def test_running_model_produces_consistent_result():
    parse_coverage_input(
        example_parameters.coverage_file_name,
        example_parameters.coverage_text_file_storage_name,
    )

    results_with_seed1 = returnYearlyPrevalenceEstimate(3.0, 0.3, seed=1, fixed_parameters=example_parameters)
    print(results_with_seed1['draw_1'])
    expected_prevalance = [0.1, 0.2, 0.3, 0.3, 0.2, 0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    pdt.assert_series_equal(results_with_seed1['draw_1'], pd.Series(expected_prevalance, name="draw_1"))

def test_running_parallel_produces_results():
    results = run_model_with_parameters(
        seeds=[1, 2, 3, 4], 
        parameters=[(3.0, 0.3), (3.0, 0.3), (3.0, 0.3), (3.0, 0.3)], 
        fixed_parameters=example_parameters,
        year_indices=[23],
        num_parallel_jobs=2)
    print(results)
    npt.assert_array_equal(results, [[0. ],[0.5],[0.4],[0.8]])

def test_running_model_with_different_seed_gives_different_result():
    parse_coverage_input(
        example_parameters.coverage_file_name,
        example_parameters.coverage_text_file_storage_name,
    )

    results_with_seed1 = returnYearlyPrevalenceEstimate(3.0, 0.3, seed=1, fixed_parameters=example_parameters)
    results_with_seed2 = returnYearlyPrevalenceEstimate(3.0, 0.3, seed=2, fixed_parameters=example_parameters)
    with pytest.raises(AssertionError):
        pdt.assert_frame_equal(results_with_seed1, results_with_seed2)

def test_extract_data():
    example_results = pd.DataFrame({"Time": [0.0, 1.0, 2.0], "draw_1": [0.1, 0.2, 0.3]})
    outcome = extract_relevant_results(example_results, [2.0])
    pdt.assert_series_equal(
        outcome, pd.Series(data=[0.3], index=[2.0], name="Prevalence")
    )


def test_extract_multiple_years():
    example_results = pd.DataFrame({"Time": [0.0, 1.0, 2.0], "draw_1": [0.1, 0.2, 0.3]})
    outcome = extract_relevant_results(example_results, [0.0, 2.0])
    pdt.assert_series_equal(
        outcome, pd.Series(data=[0.1, 0.3], index=[0.0, 2.0], name="Prevalence")
    )


def test_extract_missing_year():
    example_results = pd.DataFrame({"Time": [0.0, 1.0, 2.0], "draw_1": [0.1, 0.2, 0.3]})
    with pytest.raises(ValueError):
        outcome = extract_relevant_results(example_results, [3.0])

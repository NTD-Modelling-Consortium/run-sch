from dataclasses import dataclass
import random
import os
import time
from typing import Literal
import pandas as pd
import sch_simulation
import numpy as np
from sch_simulation.helsim_FUNC_KK.configuration import setupSD
from sch_simulation.helsim_FUNC_KK.file_parsing import (
    parse_coverage_input,
    parse_vector_control_input,
    readCoverageFile,
)

from joblib import Parallel, delayed

import sch_simulation.helsim_RUN_KK

import sch_simulation.helsim_FUNC_KK.results_processing as results_processing

@dataclass(eq=True, frozen=True)
class FixedParameters:
    number_hosts: int
    """Number of people to model

    The higher the value of N, the more consistent the results will be
    though the longer the simulation will take."""

    coverage_file_name: str
    demography_name: str
    """Which demography to use"""

    survey_type: Literal["KK1", "KK2", "POC-CCA", "PCR"]
    """Which suvery type to use
    
    Must be one KK1, KK2, POC-CCA or PCR"""

    coverage_text_file_storage_name: str
    """File name to store coverage information in"""

    parameter_file_name: str
    """Standard parameter file path (in sch_simulation/data folder)"""

    min_multiplier: int
    """Used to speed up running the model
    
    A higher number will result in faster but less accurate simulation."""


def returnYearlyPrevalenceEstimate(R0, k, seed, fixed_parameters: FixedParameters):
    np.random.seed(seed)
    # initialize the parameters
    params = sch_simulation.helsim_RUN_KK.loadParameters(
        fixed_parameters.parameter_file_name, fixed_parameters.demography_name
    )
    # add coverage data to parameters file
    params = readCoverageFile(fixed_parameters.coverage_text_file_storage_name, params)
    # add vector control data to parameters
    params = parse_vector_control_input(fixed_parameters.coverage_file_name, params)

    # update R0, k and population size as desired
    params.R0 = R0
    params.k = k
    params.N = fixed_parameters.number_hosts

    # setup the initial simData

    simData = setupSD(params)

    # run a single realization
    results, SD = sch_simulation.helsim_RUN_KK.doRealizationSurveyCoveragePickle(
        params,
        fixed_parameters.survey_type,
        simData,
        fixed_parameters.min_multiplier,
    )

    # process the output
    output = results_processing.extractHostData([results])

    # do a single simulation
    numReps = 1
    PrevalenceEstimate = results_processing.getPrevalenceWholePop(
        output, params, numReps, params.Unfertilized, fixed_parameters.survey_type, 1
    )
    return PrevalenceEstimate


def extract_relevant_results(
    results: pd.DataFrame, relevant_years: list[float]
) -> float:

    relevant_rows = results["Time"].isin(relevant_years)
    prevalence_for_relevant_years = pd.Series(
        data=results[relevant_rows][results_processing.OUTPUT_COLUMN_NAME],
        index=relevant_years,
        name="Prevalence",
    )
    if relevant_rows.sum() < len(relevant_years):
        raise ValueError(
            f"Missing data for requested years: \n{prevalence_for_relevant_years}"
        )

    return prevalence_for_relevant_years

def run_and_extract_results(parameter_set, seed, fixed_parameters, year_indices, include_output=False):
    R0 = parameter_set[0]
    k = parameter_set[1]

    start_time = time.time()

    results = returnYearlyPrevalenceEstimate(R0, k, seed, fixed_parameters)

    end_time = time.time()
    if include_output:
        print(f'Run R0: {R0}, k: {k} took {end_time-start_time:10.2f} seconds')

    return extract_relevant_results(results, year_indices)

def run_model_with_parameters(
    seeds, parameters, fixed_parameters: FixedParameters, year_indices: list[int],
    num_parallel_jobs = -2 # default to all but one process to keep computers responsive
):
    if len(seeds) != len(parameters):
        raise ValueError(
            f"Must have same number of seeds as parameters {len(seeds)} != {len(parameters)}"
        )
    
    parse_coverage_input(
        fixed_parameters.coverage_file_name,
        fixed_parameters.coverage_text_file_storage_name,
    )
    
    print(f'Running {len(seeds)} simulations across {num_parallel_jobs} cores')

    num_runs = len(seeds)

    print_timing_info_every_n_times = 10

    final_prevalence_for_each_run = Parallel(n_jobs=num_parallel_jobs)(delayed(run_and_extract_results)
        (parameter_set, seed, fixed_parameters, year_indices, include_output=index % print_timing_info_every_n_times == 0) 
        for index, (seed, parameter_set) in enumerate(zip(seeds, parameters)))

    results_np_array = np.array(final_prevalence_for_each_run).reshape(
        num_runs, len(year_indices)
    )

    os.remove(fixed_parameters.coverage_text_file_storage_name)

    return results_np_array

import argparse
import copy
import os
import pickle
import time
from pathlib import Path

from joblib import Parallel, delayed
import numpy as np
import pandas as pd

from sch_simulation.helsim_RUN_KK import (
    doRealizationSurveyCoveragePickle,
    loadParameters,
)
from sch_simulation.helsim_FUNC_KK import (
    configuration,
    file_parsing,
    results_processing,
    utils,
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run STH nearterm projections for IUs in a batch ID",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run projections for batch 10, ascaris species
  python sth_projections_per_IU.py --id=10 --species=ascaris
  
  # Run projections for batch 10, trichuris species, with 8 cores
  python sth_projections_per_IU.py --id=10 --species=trichuris --num-cores=8
"""
    )
    
    parser.add_argument(
        "--id",
        type=int,
        required=False,
        help="Batch/Task ID - processes all IUs with this TaskID. If not specified, processes all batches."
    )
    
    parser.add_argument(
        "--species",
        type=str,
        choices=["ascaris", "hookworm", "trichuris"],
        default=None,
        required=False,
        help="STH species to run projections for. If not specified, processes all species."
    )
    
    parser.add_argument(
        "--num-cores",
        type=int,
        default=10,
        help="Number of cores to use (default: 10)"
    )
    
    return parser.parse_args()


start = time.time()

"""
    extra functions to output the prevalence data
"""


def adjustR0AndKParams(parameters, R0, k):
    parameters.R0 = R0
    parameters.k = k

    # configure the parameters
    parameters = configuration.configure(parameters)
    parameters.psi = utils.getPsi(parameters)
    parameters.equiData = configuration.getEquilibrium(parameters)
    return parameters


def generateSimData(seed, params, R0, k):
    np.random.seed(seed)
    # adjust parameters here so that we will use the correct set of parameters
    # to generate the new population for each simulation
    params = adjustR0AndKParams(params, R0, k)
    # return the starting simData
    return configuration.setupSD(params)


def constructNTDMCResults(params, res, startYear):
    """
    This function will return a data frame forNTDMC containing appropriate data

    Parameters
    ----------
    params: Parameters
        dataclass containing the parameter names and values;
    res:
        the results from the simulation along with the state of the population
    Returns
    -------
     Data frame for NTDMC which contains simpler
    data, to do with prevalence of SAC and the whole population.
    """

    # unpack results from the simulation here. These will be called to extract data from later
    results = [item[0] for item in res]
    for i in range(len(results)):
        singleSimResult = results[i]

        _, _, _, dfSAC = returnNTDMCOutputs(
            params, singleSimResult, [5, 15], "SAC", startYear
        )
        _, _, _, dfAll = returnNTDMCOutputs(
            params, singleSimResult, [0, 100], "Whole Population", startYear
        )

        if i == 0:
            NTDMC = pd.concat([dfSAC, dfAll], ignore_index=True)
            NTDMC["draw_0"] = NTDMC["draw_1"].values
            NTDMC = NTDMC.drop("draw_1", axis=1)
        else:
            colname = "draw_" + str(i)

            newColsNTDMC = pd.concat([dfSAC, dfAll], ignore_index=True)
            NTDMC[colname] = newColsNTDMC["draw_1"].values

    return NTDMC


def returnNTDMCOutputs(
    params,
    results,
    ageBand,
    PopType,
    startYear,
    prevThreshold=0.02,
    surveyType="KK2",
    numReps=1,
    nSamples=2,
    sampleSize=100,
):
    output = results_processing.extractHostData([results])
    prevalence, _, medium_prevalence, heavy_prevalence, _ = (
        results_processing.getBurdens(
            output,
            params,
            numReps,
            ageBand,
            params.Unfertilized,
            surveyType,
            nSamples,
            sampleSize,
        )
    )
    allTimes = output[0].timePoints + startYear
    if PopType == "SAC":
        prevBelowThreshold = (medium_prevalence + heavy_prevalence) < prevThreshold
    newrows = pd.DataFrame(
        {
            "year_id": allTimes,
            "age_start": np.repeat(ageBand[0], len(allTimes)),
            "age_end": np.repeat(ageBand[1], len(allTimes)),
            "intensity": np.repeat("None", len(allTimes)),
            "species": np.repeat(params.species, len(allTimes)),
            "measure": np.repeat("Prevalence " + PopType, len(allTimes)),
            "draw_1": prevalence,
        }
    )
    df1 = newrows
    newrows = pd.DataFrame(
        {
            "year_id": allTimes,
            "age_start": np.repeat(ageBand[0], len(allTimes)),
            "age_end": np.repeat(ageBand[1], len(allTimes)),
            "intensity": np.repeat("None", len(allTimes)),
            "species": np.repeat(params.species, len(allTimes)),
            "measure": np.repeat("Medium + Heavy Prevalence " + PopType, len(allTimes)),
            "draw_1": medium_prevalence + heavy_prevalence,
        }
    )
    df1 = pd.concat([df1, newrows], ignore_index=True)
    if PopType == "SAC":
        newrows = pd.DataFrame(
            {
                "year_id": allTimes,
                "age_start": np.repeat(ageBand[0], len(allTimes)),
                "age_end": np.repeat(ageBand[1], len(allTimes)),
                "intensity": np.repeat("None", len(allTimes)),
                "species": np.repeat(params.species, len(allTimes)),
                "measure": np.repeat("Below  EPHP threshold", len(allTimes)),
                "draw_1": prevBelowThreshold,
            }
        )
        df1 = pd.concat([df1, newrows], ignore_index=True)
    return prevalence, medium_prevalence, heavy_prevalence, df1


############################################################################################################
############################################################################################################
############################################################################################################

def run_projections_for_iu(iu, country, species, species_prefix, args):
    """Run projections for a single IU."""
    print(f"\nProcessing IU {iu} in {country} for species {species}")
    
    # Get environment paths
    path_to_inputs = Path(os.getenv("PATH_TO_PROJECTIONS_INPUTS", "."))
    path_to_artefacts = Path(os.getenv("PATH_TO_PROJECTIONS_ARTEFACTS", "."))
    path_to_fitting_prep_artefacts = Path(os.getenv("PATH_TO_FITTING_PREP_ARTEFACTS", "."))
    path_to_projections_prep_artefacts = Path(os.getenv("PATH_TO_PROJECTIONS_PREP_ARTEFACTS", "."))
    path_to_model = Path(os.getenv("STH_SCH_MODEL_DIR", "/ntdmc/sth-sch-amis-integration/model/ntd-model-sch"))
    
    # Coverage file path (from projections-prep artefacts)
    if species == "trichuris":
        coverage_file_name = f"endgame_inputs/STH/InputMDA_MTP_projections_trichuris_{iu}.xlsx"
    else:
        coverage_file_name = f"endgame_inputs/STH/InputMDA_MTP_projections_{iu}.xlsx"
    
    coverage_file_path = path_to_fitting_prep_artefacts / coverage_file_name
    
    if not os.path.exists(coverage_file_path):
        print(f"Warning: Coverage file not found: {coverage_file_path}")
        print(f"Skipping IU {iu}")
        return False
    
    # File paths
    demog_name = "UgandaRural"
    coverage_text_file_storage_name = path_to_projections_prep_artefacts / f"Man_MDA_vacc/Man_MDA_vacc_{species}_{iu}.txt"
    # Parameter file is in the model package - loadParameters expects just the filename relative to data/
    param_file_name = f"STH_params/{species}_params_projections.txt"
    
    # Parameters file path
    rk_file_path = path_to_projections_prep_artefacts / f"InputPars_MTP_{species}/InputPars_MTP_{iu}.csv"
    
    if not os.path.exists(rk_file_path):
        print(f"Warning: Parameters file not found: {rk_file_path}")
        print(f"Skipping IU {iu}")
        return False
    
    start_year = 1985
    num_sims = 200
    
    # Read in parameter and coverage files
    _ = file_parsing.parse_coverage_input(coverage_file_path, str(coverage_text_file_storage_name))
    # Initialize the parameters
    params = loadParameters(param_file_name, demog_name)
    # Add coverage data to parameters file
    params = file_parsing.readCoverageFile(str(coverage_text_file_storage_name), params)
    # Add vector control data to parameters
    params = file_parsing.parse_vector_control_input(coverage_file_path, params)
    
    # Read in fitted parameters for IU
    simparams = pd.read_csv(rk_file_path)
    simparams.columns = [s.replace(" ", "") for s in simparams.columns]
    
    # Define the lists of random seeds, R0 and k
    seed = simparams.iloc[:, 1].tolist()
    R0 = simparams.iloc[:, 2].tolist()
    k1 = simparams.iloc[:, 3].tolist()
    
    # Setup the output times to be every year
    params.outTimings = np.arange(0, params.maxTime, 1)
    
    # Set the survey type to Kato Katz with duplicate slide
    survey_type = "KK2"
    
    params_to_alter = copy.deepcopy(params)
    
    # Run the simulations
    print(f"Running {num_sims} simulations for IU {iu}...")
    res = Parallel(n_jobs=args.num_cores)(
        delayed(doRealizationSurveyCoveragePickle)(
            adjustR0AndKParams(params_to_alter, R0[i], k1[i]),
            survey_type,
            generateSimData(seed[i], params_to_alter, R0[i], k1[i]),
        )
        for i in range(num_sims)
    )
    
    simData = [item[1] for item in res]
    
    # Create output directory structure
    output_dir = path_to_artefacts / "projections" / species / country / f"{country}{str(iu).zfill(5)}"
    os.makedirs(str(output_dir), exist_ok=True)

    # Output pickle file, we want outputs like <ascaris-folder>/AGO/AGO02049/Asc_AGO02049.p
    pickle_file_path = output_dir / f"{species_prefix}{country}{str(iu).zfill(5)}.p"
    print(f"Saving pickle file: {pickle_file_path}")
    pickle.dump(simData, open(str(pickle_file_path), "wb"))
    
    # Output prevalence dataset
    NTDMC = constructNTDMCResults(params, res, start_year)
    prev_dataset_file_path = output_dir / f"PrevDataset_{species_prefix}{country}{str(iu).zfill(5)}.csv"
    print(f"Saving prevalence dataset: {prev_dataset_file_path}")
    NTDMC.to_csv(str(prev_dataset_file_path), index=False)
    
    print(f"✓ Finished projections for {species} in IU {iu}")
    return True


def main():
    """Main function to run STH projections for a batch of IUs."""
    args = parse_arguments()
    
    # Get task ID from argument - if not provided, process all batches
    process_single_batch = args.id is not None
    if process_single_batch:
        task_ids_to_process = [args.id]
    else:
        # Will determine all task IDs from the lookup table
        task_ids_to_process = None
    
    species = args.species
    
    # Species prefix mapping
    species_prefix_map = {
        "ascaris": "Asc_",
        "hookworm": "Hook_",
        "trichuris": "Tri_"
    }
    species_prefix = species_prefix_map[species]
    
    print(f"Using {args.num_cores} cores for parallel processing")
    
    # Get path to lookup table
    path_to_fitting_prep = os.getenv("PATH_TO_FITTING_PREP_ARTEFACTS", ".")
    lookup_file_path = os.path.join(path_to_fitting_prep, "Maps/table_iu_idx_STH.csv")
    
    if not os.path.exists(lookup_file_path):
        raise FileNotFoundError(f"Lookup table not found: {lookup_file_path}")
    
    # Read lookup table
    df_iu_country = pd.read_csv(lookup_file_path)
    
    # Determine which task IDs to process
    if task_ids_to_process is None:
        # Process all batches
        task_ids_to_process = sorted(df_iu_country["TaskID"].unique())
        print(f"Running STH projections for ALL TaskIDs ({len(task_ids_to_process)} batches), species: {species}")
        print(f"TaskIDs to process: {task_ids_to_process[:10]}{'...' if len(task_ids_to_process) > 10 else ''}")
    else:
        # Process specific batch
        task_id = task_ids_to_process[0]
        batch_ius = df_iu_country[df_iu_country["TaskID"] == task_id]
        if batch_ius.empty:
            raise ValueError(f"No IUs found for TaskID {task_id}")
        print(f"Running STH projections for TaskID {task_id}, species: {species}")
        print(f"Found {len(batch_ius)} IUs for TaskID {task_id}")
    
    # Process each task ID
    total_successful = 0
    total_failed = 0
    
    for task_id in task_ids_to_process:
        if len(task_ids_to_process) > 1:
            print(f"\n{'='*40}")
            print(f"Processing TaskID {task_id}")
            print(f"{'='*40}")
        
        batch_ius = df_iu_country[df_iu_country["TaskID"] == task_id]
        successful_ius = []
        failed_ius = []
        
        for _, row in batch_ius.iterrows():
            iu = row["IU_CODE"]
            country = row["country"]
            
            try:
                if run_projections_for_iu(iu, country, species, species_prefix, args):
                    successful_ius.append(iu)
                else:
                    failed_ius.append(iu)
            except Exception as e:
                print(f"✗ Error processing IU {iu}: {e}")
                failed_ius.append(iu)
        
        total_successful += len(successful_ius)
        total_failed += len(failed_ius)
        
        if len(task_ids_to_process) > 1:
            print(f"TaskID {task_id}: ✓ {len(successful_ius)} successful, ✗ {len(failed_ius)} failed")
    
    # Final summary
    print(f"\n{'='*60}")
    if len(task_ids_to_process) == 1:
        print(f"PROJECTION SUMMARY for TaskID {task_ids_to_process[0]}, species: {species}")
    else:
        print(f"FINAL PROJECTION SUMMARY - species: {species}")
    print(f"{'='*60}")
    print(f"✓ Total successful IUs: {total_successful}")
    if total_failed > 0:
        print(f"✗ Total failed IUs: {total_failed}")
    print(f"Total processing time: {time.time() - start:.2f} seconds")
    
    if total_failed > 0:
        raise RuntimeError(f"Projections failed for {total_failed} IUs")
    
    print(f"\n✓ All projections completed successfully!")


if __name__ == "__main__":
    main()

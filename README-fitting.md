Scripts used for STH fitting and near term projections 
================

Notes about SCH: 
- Warning: I adjusted the code to match this repo but it hasn't been tested. 
- For mansoni: when last ran I did the full set of results (low and high burden) for all IUs and then produced a CSV that said which version to use for each IU (the one with the best ESS) for the projections to 2040. But when redoing should change this to use the model evidence (wasn't available at the time of last runs).
- SCH also wasn't fitted with importation so the parameter files in `ntd-model-sch/sch_simulation/data/SCH_params` need to be updated if we want to use importation


### Installation

- Install virtual environment. Following instructions are for the Oxford cluster

  ```
  module purge
  module load Python/3.11.3-GCCcore-12.3.0
  python3 -m venv .venv
  source .venv/bin/activate
  ```
  
- At the time of running the scripts should be run using the branch `updateImportation` found of the `ntd-model-sch` repo <https://github.com/NTD-Modelling-Consortium/ntd-model-sch/>.

  ```
  git clone git@github.com:NTD-Modelling-Consortium/ntd-model-sch.git
  cd ntd-model-sch
  git fetch
  git switch updateImportation
  pip install .
  ```
  
- In order for the model to talk to AMIS, use the `run-sch` repo <https://github.com/NTD-Modelling-Consortium/run-sch/>. Deceptively named, this can run both STH and SCH since they use the same model.

  ```
  cd ..
  git clone git@github.com:NTD-Modelling-Consortium/run-sch.git
  cd run-sch
  pip install .
  ```

- Make sure required R libraries are installed. The `run-sch` repo has `renv.lock` files for different versions of R. But older versions of these files may be out of date and using the non-CRAN version of `AMISforInfectiousDiseases`. I don't think I actually used these, I just installed the required libraries directly (I think it's just dplyr, AMISforInfectiousDiseases, mvtnorm). For prep of the maps/histories we also require tidyr, writexl, readxl, pracma. For the Oxford cluster the version used was:

  ```
  module load R/4.3.2-gfbf-2023a
  ```

- On HPC clusters, scripts that take long to run must be run through Slurm. Tables below show shell scripts because of this.

### Preparing histories and maps for the fitting 

| R script  (see `Maps-STH/` and `Maps-SCH/` directories)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| prepare_histories.R                                         | run_hist_maps.sh              |
| prepare_histories_trichuris.R                               | run_hist_maps.sh              |
| prepare_maps_all_species.R                                  | run_hist_maps.sh              |
| prepare_histories_projections.R                             | run_prep_hist_proj.sh         |
| prepare_histories_trichuris_projections.R                   | run_prep_hist_proj.sh         |
| prepare_histories_and_maps_mansoni.R                   | runPrep_fitting.sh         |
| prepare_histories_and_maps_haematobium.R                   | runPrep_fitting.sh         |
| prepare_histories_projections_sch.R                   | runPrep_projections.sh         |

- If the data/batch allocations are changing for SCH need to manually change `id_no_mda` in `prepare_histories_projections_sch.R` to reflect batches with no MDA
- Similarly for trichuris, if data/batch allocations change then need to change `original_last_batch_ID` in `prepare_histories_trichuris_projections`
- Not required for changes to ascaris/hookworm because its due to the batches being reassigned at some point in the respective preparation scripts

### Running the fitting

| R script  (see `run-sch/` directory)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| sth_fitting.R                                               | runFit_asca.sh, runFit_hook.sh, runFit_tric.sh  |            |
| sth_fitting_sigma0.025.R                                    | runFit_asca_sigma0.025.sh, runFit_hook_sigma0.025.sh, runFit_tric_sigma0.025.sh             |
| sch_fitting.R                                               | runFit_asca.sh, runFit_hook.sh, runFit_tric.sh  |            |
| sch_fitting_sigma0.025.R                                    | runFit_asca_sigma0.025.sh, runFit_hook_sigma0.025.sh, runFit_tric_sigma0.025.sh             |
| find_lowESS_ids.R                                           | runFindLowESS.sh         |


- Before running these, we have to manually set the corresponding `species` in 
`run-sch/sth_fitting.R`, `run-sch/sth_fitting_sigma0.025.R` and `find_lowESS_ids.R`
- You will also need to specify `failed_ids` (i.e. batches that failed during sigma=0.0025 runs) in `find_lowESS_ids.R` (use something like `grep -i "error" <log file names>` to find failed batches)
- `sth_fitting_sigma0.025.R` should be run for batches that: are in `failed_ids` or the output printed to console (or log files if using the shell script) from `find_lowESS_ids.R` (by species)
- Note: the failed batches probably contain some IUs that are actually able to be fit, but the fitting failed after these dropped out of the active set, or the fitting timed out, and I didn't have time to go back and refit these

### Prepare for near term projections

These are for until the end of 2025, which means 2026.0 in the continuous scale.

| R script  (see `post_AMIS_analysis/` directory)             | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| preprocess_for_projections.R                                | runPreprocessing.sh           |

<br/>
**preprocess_for_projections.R**: creates the 200 parameter vectors (simulated from the fitted 
models) used in projections. Reorganises the files with the 200 samples used in projections, 
so that they are organised in the expected file hierarchy in the cloud.

- Before running `runPreprocessing.sh`, we have to manually choose `species` in `post_AMIS_analysis/preprocessing_for_projections.R` 
- You also need to specify in `post_AMIS_analysis/preprocessing_for_projections.R` the batches that failed **when sigma=0.0025 and sigma=0.025** (`failed_ids` and `failed_ids_sigma0.025` respectively)
- Note this relies on outputs from `Maps/prepare_histories_projections*` files

### Plots for the model fits

| R script  (see `post_AMIS_analysis/` directory)             | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| out_all_countries_STH.R                                         | NA      |
| IUsWithInsufficientESS.R                                    | NA                            |     


<br/>
- **out_all_countries_STH.R**:  saves plots in and summary maps for STH
- **IUsWithInsufficientESS.R**:  finds IUs that have ESS < 200 (after also trying higher sigma=0.025). 
- In both of these files you also need to specify the batches that failed **when sigma=0.0025 and sigma=0.025** (`failed_ids` and `failed_ids_sigma0.025` respectively)

### Running the near term projections

| Python script  (see `run-sch/sth_amis/` directory)      | Corresponding shell script    |
|:--------------------------------------------------------|:------------------------------|
| sth_projections_per_IU.py                               | runProj_STH.sh                |
| sch_projections_per_IU.py                               | runProj_SCH.sh                |

- There can be a maximum number of tasks that can be submitted at a time on HPC clusters. 

- Before running **runProj_STH.sh**, we have to manually choose `species` in `run-sch/sth_amis/sth_projections_per_IU.py`. Similiarly for **runProj_SCH.sh**

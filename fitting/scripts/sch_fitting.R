library(dplyr)
library(AMISforInfectiousDiseases)
library(optparse)

# Get paths from environment variables
kPathToFittingPrepArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToMaps <- file.path(kPathToFittingPrepArtefacts, "Maps")
kPathToArtefacts <- Sys.getenv("PATH_TO_FITTING_ARTEFACTS")
kPathToTrajectories <- file.path(kPathToArtefacts, "trajectories")
kPathToAmisOutput <- file.path(kPathToArtefacts, "AMIS_output")
kPathToInfections <- file.path(kPathToArtefacts, "infections")
kPathToFittingScripts <- Sys.getenv("PATH_TO_FITTING_SCRIPTS")

# Create output directories if they don't exist
if (!dir.exists(kPathToTrajectories)) {
  dir.create(kPathToTrajectories, recursive = TRUE)
}
if (!dir.exists(kPathToAmisOutput)) {
  dir.create(kPathToAmisOutput, recursive = TRUE)
}
if (!dir.exists(kPathToInfections)) {
  dir.create(kPathToInfections, recursive = TRUE)
}

# Command line arguments
option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    help = "Single batch ID to process. If not provided, will check the 'SLURM_ARRAY_TASK_ID' environment variable."
  ),
  make_option(c("-s", "--species"),
    type = "character",
    help = "Species to process: haematobium, mansoni_low_burden, or mansoni_high_burden",
    metavar = "SPECIES"
  ),
  make_option(c("--amis-sigma"),
    type = "double",
    default = 0.0025,
    help = "AMIS sigma parameter (default: 0.0025)"
  ),
  make_option(c("--amis-n-samples"),
    type = "integer",
    default = 500,
    help = "Number of AMIS samples (default: 500)"
  ),
  make_option(c("--amis-target-ess"),
    type = "integer",
    default = 500,
    help = "Target ESS parameter for AMIS (default: 500)"
  ),
  make_option(c("--num-cores"),
    type = "integer",
    default = NULL,
    help = "Number of cores to use for parallel processing (default: all available)"
  ),
  make_option(c("--amis-n-iters"),
    type = "integer",
    default = 50,
    help = "Maximum number of AMIS iterations (default: 50)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Get batch ID
id <- if (!is.null(opts$id)) {
  opts$id
} else {
  env_var_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if (!is.na(env_var_id)) {
    env_var_id
  } else {
    stop("Batch ID is not provided in the command line arguments and SLURM_ARRAY_TASK_ID environment variable is not defined")
  }
}

# Get species
if (is.null(opts$species)) {
  stop("Species must be specified using --species argument (haematobium, mansoni_low_burden, or mansoni_high_burden)")
}
species <- opts$species
if (!species %in% c("haematobium", "mansoni_low_burden", "mansoni_high_burden")) {
  stop("Species must be one of: haematobium, mansoni_low_burden, mansoni_high_burden")
}

# Get number of cores
num_cores_to_use <- if (!is.null(opts$"num-cores")) {
  opts$"num-cores"
} else {
  parallel::detectCores()
}

print(paste0("Processing batch ID ", id, " for species ", species))
print(paste0("Using ", num_cores_to_use, " cores"))
print(paste0("AMIS parameters: sigma=", opts$"amis-sigma", ", n_samples=", opts$"amis-n-samples", ", target_ess=", opts$"amis-target-ess", ", max_iters=", opts$"amis-n-iters"))

# Source AMIS integration code
source(file.path(kPathToFittingScripts, "amis_integration.R"))

# Setup Python environment
# In Docker, Python packages are installed globally, no venv needed
# reticulate::use_virtualenv("../.venv", required=TRUE)
sch_simulation <- get_amis_integration_package()
reticulate::py_config()

# Determine endgame inputs directory based on species
endgame_dir <- if (species == "haematobium") {
  "sch-haematobium"
} else {
  "sch-mansoni"
}

# Fixed parameters for the simulation
fixed_parameters <- sch_simulation$FixedParameters(
  # the higher the value of N, the more consistent the results will be
  # though the longer the simulation will take
  number_hosts = 500L,
  # MDA coverage file
  coverage_file_name = file.path(kPathToFittingPrepArtefacts, "endgame_inputs", endgame_dir, paste0("InputMDA_MTP_", id, ".xlsx")),
  demography_name = "UgandaRural",
  # set the survey type to Kato Katz with duplicate slide
  survey_type = "KK2",
  parameter_file_name = paste0("SCH_params/", species, "_params.txt"),
  coverage_text_file_storage_name = paste0("Man_MDA_vacc_", species, "_", id, ".txt"),
  # the following number dictates the number of events (e.g. worm deaths)
  # we allow to happen before updating other parts of the model
  # the higher this number the faster the simulation
  # (though there is a check so that there can't be too many events at once)
  # the higher the number the greater the potential for
  # errors in the model accruing.
  # 5 is a reasonable level of compromise for speed and errors, but using
  # a lower value such as 3 is also quite good
  min_multiplier = 5L
)

# Clarifying year_indices
# - The Python code can simulate prevalences at any point in continuous time.
# - In the Python code, year_indices 0, 1, ..., 33 refer to years 1985.0, 1986.0, ..., 2018.0 (beginning of each year).
# - Map samples are in discrete time, and we assume that samples for 2002 refer to the end of the year 2002 (approximately, 2002.99999 ~ 2003.0).
# - Therefore, if we have map samples for the (end of) years 2002, 2013, 2022, these should be compared to the simulations at times 2003.0, 2014.0, 2023.0 if we use the discrete-valued vector "year_indices".
# - Thus, we need to pass year_indices <- c(18L,29L,38L) to the Python code
# - The Python model will return a matrix with only 3 columns called 18L,29L,38L.
year_indices <- c(18L, 29L, 38L)

# Load prevalence map based on species type
if (species %in% c("mansoni_low_burden", "mansoni_high_burden")) {
  load(file.path(kPathToMaps, "mansoni_maps.rds"))
  prevalence_map <- mansoni_maps
} else {
  load(file.path(kPathToMaps, paste0(species, "_maps.rds")))
  prevalence_map <- get(paste0(species, "_maps"))
}

# Filter for TaskID
prevalence_map <- lapply(1:length(prevalence_map), function(t) {
  output <- list(data = as.matrix(prevalence_map[[t]]$data %>%
                                   filter(TaskID == id) %>%
                                   select(-c(IU_ID, TaskID))))
  rownames(output$data) <- prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID == id]
  return(output)
})

# Load prior
source(file.path(kPathToFittingScripts, "sch_prior.R"))
prior <- Prior

# Algorithm parameters
amis_params <- default_amis_params()
amis_params$max_iters <- opts$"amis-n-iters"
amis_params$n_samples <- opts$"amis-n-samples"
amis_params$target_ess <- opts$"amis-target-ess"
amis_params$sigma <- opts$"amis-sigma"
amis_params$boundaries <- c(0, 1)

# Shell to save trajectories
trajectories <- c() # save simulated trajectories as code is running
path_trajectories_file <- file.path(kPathToTrajectories, paste0("trajectories_", id, "_", species, ".Rdata"))
save(trajectories, file = path_trajectories_file)

# Run AMIS
st <- Sys.time()
amis_output <- AMISforInfectiousDiseases::amis(
  prevalence_map,
  build_transmission_model(
    prevalence_map, fixed_parameters, year_indices, num_cores_to_use, NULL,
    path_trajectories_file
  ),
  prior,
  amis_params,
  seed = id
)
en <- Sys.time()
dur_amis <- as.numeric(difftime(en, st, units = "mins"))

# Save AMIS output with sigma suffix if not default
output_suffix <- ifelse(opts$"amis-sigma" == 0.0025, "", paste0("_sigma", opts$"amis-sigma"))
save(amis_output, file = file.path(kPathToAmisOutput, paste0(species, "_amis_output", id, output_suffix, ".Rdata")))

print(amis_output)
summary(amis_output)
cat("--------------------- \n")
print(paste0("AMIS run time: ", round(dur_amis, digits = 2), " minutes"))

# Save infections
infections <- rowMeans(t(prevalence_map[[3]]$data))
save(infections, file = file.path(kPathToInfections, paste0("infections", id, "_", species, ".Rdata")))

# Save summary
ess <- amis_output$ess
n_success <- length(which(ess >= amis_params[["target_ess"]]))
failures <- which(ess < amis_params[["target_ess"]])
n_failure <- length(failures)

# Save ESS failures
ess_file <- file.path(kPathToArtefacts, paste0("ESS_NOT_REACHED_", species, output_suffix, ".txt"))
if (n_failure > 0) {
  cat(paste(failures, id, ess[failures]), file = ess_file, sep = "\n", append = TRUE)
}

# Save summary
summary_file <- file.path(kPathToArtefacts, paste0("summary_", species, output_suffix, ".csv"))
if (!file.exists(summary_file)) {
  cat("ID,n_failure,n_success,n_sim,min_ess,duration_amis,duration_subsampling\n", file = summary_file)
}
cat(id, n_failure, n_success, length(amis_output$seeds), min(ess), dur_amis, NA, "\n",
  sep = ",", file = summary_file, append = TRUE
)

cat("Fitting completed successfully!\n")

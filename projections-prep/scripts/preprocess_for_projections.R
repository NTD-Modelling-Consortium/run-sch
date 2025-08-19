library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")
library("optparse")

# Define command line options
option_list <- list(
  make_option(c("-s", "--species"),
    type = "character",
    default = "ascaris",
    help = "Species to process (ascaris, hookworm, trichuris, haematobium, mansoni_low_burden, mansoni_high_burden) [default=%default]"
  ),
  make_option(c("-i", "--id"),
    type = "integer",
    default = NULL,
    help = "Single batch ID to process. If not provided, will process all batches"
  ),
  make_option(c("-f", "--failed-ids"),
    type = "character",
    default = NULL,
    help = "Comma-separated list of failed batch IDs that used sigma=0.025"
  ),
  make_option(c("--seed"),
    type = "integer",
    default = 100,
    help = "Random seed [default=%default]"
  ),
  make_option(c("--amis-sigma"),
    type = "double",
    default = 0.0025,
    help = "AMIS sigma parameter (default: 0.0025)"
  ),
  make_option(c("--ess-threshold"),
    type = "integer",
    default = 200,
    help = "ESS threshold parameter (default: 200)"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Set random seed
set.seed(opts$seed)

# Get species from command line
species <- opts$species

# Validate species choice
valid_species <- c("ascaris", "hookworm", "trichuris", "haematobium", "mansoni_low_burden", "mansoni_high_burden")
if (!species %in% valid_species) {
  stop(paste("Invalid species:", species, ". Valid choices:", paste(valid_species, collapse=", ")))
}

# Parse failed IDs from command line or use defaults
failed_ids <- c()
if (!is.null(opts$"failed-ids")) {
  failed_ids <- as.numeric(strsplit(opts$"failed-ids", ",")[[1]])
}

# Set up paths from environment variables
kPathToFittingPrepArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToMapsArtefacts <- file.path(kPathToFittingPrepArtefacts, "Maps")
kPathToFittingArtefacts <- Sys.getenv("PATH_TO_FITTING_ARTEFACTS")
kPathToProjectionsPrepArtefacts <- Sys.getenv("PATH_TO_PROJECTIONS_PREP_ARTEFACTS")

if (kPathToFittingPrepArtefacts == "" || kPathToFittingArtefacts == "" || kPathToProjectionsPrepArtefacts == "") {
  stop("Environment variables PATH_TO_FITTING_PREP_ARTEFACTS, PATH_TO_FITTING_ARTEFACTS and PATH_TO_PROJECTIONS_PREP_ARTEFACTS must be set")
}

# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
if(species == "trichuris"){
  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_trichuris.rds"))
}else if (species %in% c("ascaris","hookworm")){
  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_sth.rds"))  
}else if (species == "haematobium"){
  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_haema.rds"))
} else {
  # For mansoni, use the species-specific filename with variant
  lookup_file <- paste0("iu_task_lookup_", species, ".rds")
  load(file.path(kPathToMapsArtefacts, lookup_file))
}
 
num_batches <- max(iu_task_lookup$TaskID)

cat(paste0("species: ",species, " \n"))
cat(paste0("num_batches: ",num_batches, " \n"))

# Determine which batch IDs to process
if (!is.null(opts$id)) {
  # Single batch mode
  ids_sample_pars <- opts$id
  cat(paste0("Processing single batch ID: ", opts$id, "\n"))
} else {
  # All batches mode - skip the failed batches
  ids_sample_pars = setdiff(1:num_batches, failed_ids)
  cat(paste0("Processing all batches except failed IDs: ", paste(failed_ids, collapse=", "), "\n"))
}


# Create required directories using environment variable paths
if (!dir.exists(file.path(kPathToProjectionsPrepArtefacts, "Man_MDA_vacc"))) {
  dir.create(file.path(kPathToProjectionsPrepArtefacts, "Man_MDA_vacc"), recursive = TRUE)
}

# Directory to save sampled parameters
InputPars_MTP_path_species <- file.path(kPathToProjectionsPrepArtefacts, paste0("InputPars_MTP_", species))
if (!dir.exists(InputPars_MTP_path_species)) {
  dir.create(InputPars_MTP_path_species, recursive = TRUE)
}

# sample parameters and save draws
sampled_params_all = c()
for(id in ids_sample_pars){

  #### Load AMIS output (this is just to get IU names and initial ESS)
  amis_output_path <- file.path(kPathToFittingArtefacts, "AMIS_output")
  if(!id %in% failed_ids){
    load(file.path(amis_output_path, paste0(species, "_amis_output", id, ".Rdata"))) # loads amis_output
  } else {
    load(file.path(amis_output_path, paste0(species, "_amis_output", id, "_sigma", format(opts$"amis-sigma", scientific = FALSE), ".Rdata"))) # loads amis_output
  }

  ess = amis_output$ess
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)

  if(!id %in% failed_ids){
    iu_names_lt_ess_threshold = iu_names[ess<opts$"ess-threshold"]
    iu_names_ge_ess_threshold = iu_names[!ess<opts$"ess-threshold"]
  } else {
    iu_names_lt_ess_threshold = iu_names # if failed when sigma=0.0025 then use sigma=opts$amis-sigma for all IUs
    iu_names_ge_ess_threshold = NULL
  }

  #### Sample draws from the posterior
  num_sub_samples_posterior <- 200
	  
  for (iu in iu_names) {
    # use sigma=amis-sigma results if ESS < ess-threshold when sigma=default(0.0025)
    # If we're working with a failed batch or the IU is below ess-threshold,
    # then we'll use the amis_output computed with a different sigma (!=default(0.0025))
    if((!id %in% failed_ids) & (iu %in% iu_names_ge_ess_threshold)){
      load(file.path(amis_output_path, paste0(species, "_amis_output", id, ".Rdata"))) # loads amis_output
    } else {
      load(file.path(amis_output_path, paste0(species, "_amis_output", id, "_sigma", format(opts$"amis-sigma", scientific = FALSE), ".Rdata"))) # loads amis_output
    }
    
    sampled_params <- sample_parameters(x = amis_output, n_samples = num_sub_samples_posterior, locations = which(iu_names==iu))

    # this puts all ius in the same folder
    file_name <- file.path(InputPars_MTP_path_species, paste0("InputPars_MTP_", iu, ".csv"))
    write.csv(sampled_params, file=file_name, row.names = F)

    sampled_params_iu = cbind(IU_ID=iu,sampled_params)
    sampled_params_all = rbind(sampled_params_all,sampled_params_iu)

  }

  if(id%%100==0){cat(paste0("id=",id, "; "))}

}

save(sampled_params_all, file = file.path(InputPars_MTP_path_species, paste0("InputPars_MTP_allIUs_", species, ".rds")))
cat(paste0("Produced samples for all IUs in InputPars_MTP_",species, "/ \n"))


# Realocate files in correct format for projections

if(species %in% c("ascaris","hookworm","trichuris")){
  df_IU_country <- read.csv(file.path(kPathToMapsArtefacts, "table_iu_idx_STH.csv")) # same for all species as just want country codes
} else if (species == "haematobium"){ 
  df_IU_country <- read.csv(file.path(kPathToMapsArtefacts, "table_iu_idx_haematobium.csv"))
} else if (species %in% c("mansoni_low_burden", "mansoni_high_burden")) {
  # Handle mansoni variants with species-specific files
  df_IU_country <- read.csv(file.path(kPathToMapsArtefacts, paste0("table_iu_idx_", species, ".csv")))
} else {
  # Fallback for generic mansoni
  df_IU_country <- read.csv(file.path(kPathToMapsArtefacts, "table_iu_idx_mansoni.csv"))
}

countries <- sort(unique(df_IU_country$country))

# Create directory structure for projections
proj <- file.path(kPathToProjectionsPrepArtefacts, "projections")
if (!dir.exists(proj)) {dir.create(proj, recursive = TRUE)}

path_species <- file.path(proj, species)
if (!dir.exists(path_species)) {dir.create(path_species, recursive = TRUE)}

for(country in countries){
  path_country <- file.path(path_species, country)
  if (!dir.exists(path_country)) {dir.create(path_country, recursive = TRUE)}
}
  
# prefix
if(species == "haematobium"){
  species_prefix <- "Haema_"
}else if(species == "mansoni_high_burden"){
  species_prefix <- "Man_High_"
}else if(species == "mansoni_low_burden"){
  species_prefix <- "Man_Low_" 
}else if(species == "ascaris"){
  species_prefix <- "Asc_"
}else if(species == "hookworm"){
  species_prefix <- "Hook_"
}else if(species == "trichuris"){
  species_prefix <- "Tri_" 
}

# realocate
for(id in ids_sample_pars){
    
  if (species %in% c("ascaris","hookworm","trichuris")){
    iu_names <- as.character(unique(sort(iu_task_lookup$IU_2021[iu_task_lookup$TaskID==id])))
  }else{
    iu_names <- as.character(unique(sort(iu_task_lookup$IU_ID[iu_task_lookup$TaskID==id])))
  }
  num_samples <- 200

  for (iu in iu_names) {
    
    wh <- which(df_IU_country$IU_CODE==iu)
    if(length(wh)!=1){stop("iu must be found exactly once in df_IU_country")}
    country <- df_IU_country[wh, "country"]
    
    iu0 <- sprintf("%05d", as.integer(iu))
    path_iu <- file.path(path_species, country, paste0(country, iu0))
    if (!dir.exists(path_iu)) {dir.create(path_iu, recursive = TRUE)}
    
    file_name_old <- file.path(InputPars_MTP_path_species, paste0("InputPars_MTP_", iu, ".csv"))
    sampled_params <- read.csv(file_name_old)
    
    # Select appropriate columns based on species type
    if(species %in% c("ascaris", "hookworm", "trichuris")) {
      # STH species use "R0" (no underscore)  
      sampled_params <- sampled_params[, c("seed", "R0", "k")]
    } else {
      # SCH species use "R_0" (with underscore)
      sampled_params <- sampled_params[, c("seed", "R_0", "k")]
    }
    colnames(sampled_params) <- c("seed", "r0", "k")
    file_name_new <- file.path(path_iu, paste0("Input_Rk_", species_prefix, country, iu0, ".csv"))
    write.csv(sampled_params, file=file_name_new, row.names = F)
    
  }
  
  if(id%%100==0){cat(paste0("id=",id, "; "))}

}

cat(paste0("Samples realocated for all IUs in InputPars_MTP_", species, "/ \n"))
#!/usr/bin/env Rscript

# Find batches with low ESS (Effective Sample Size)
# Compatible with STH/SCH AMIS Integration Pipeline Docker environment

library("optparse")
library("AMISforInfectiousDiseases")
library("dplyr")

# Command line arguments
option_list <- list(
  make_option(c("-s", "--species"), type="character", default=NULL,
              help="Species to analyze (ascaris, hookworm, trichuris, haematobium, mansoni_low_burden, mansoni_high_burden)", metavar="character"),
  make_option(c("-t", "--ess-threshold"), type="double", default=200,
              help="ESS threshold below which batches are flagged [default %default]", metavar="number"),
  make_option(c("-f", "--failed-ids"), type="character", default="",
              help="Comma-separated list of batch IDs that completely failed [default none]", metavar="character"),
  make_option(c("--sigma"), type="character", default="0.0025",
              help="Sigma value for AMIS output files [default %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, 
                          description="Analyze AMIS fitting results to identify batches with insufficient ESS")
opts <- parse_args(opt_parser)

if(is.null(opts$species)){
  print_help(opt_parser)
  stop("Species must be specified with --species")
}

species <- opts$species
ess_threshold <- opts$ess_threshold
sigma_value <- opts$sigma

# Parse failed IDs
if(opts$`failed-ids` == "") {
  failed_ids <- c()
} else {
  failed_ids <- as.numeric(strsplit(opts$`failed-ids`, ",")[[1]])
}

cat(paste0("Analyzing species: ", species, "\n"))
cat(paste0("ESS threshold: ", ess_threshold, "\n"))
cat(paste0("Sigma value: ", sigma_value, "\n"))
if(length(failed_ids) > 0) {
  cat(paste0("Failed batch IDs to exclude: ", paste(failed_ids, collapse=", "), "\n"))
}
cat("\n")

# Use environment variables for paths (Docker-compatible)
kPathToFittingArtefacts <- Sys.getenv("PATH_TO_FITTING_ARTEFACTS", ".")
kPathToMapsArtefacts <- file.path(Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS", "."), "Maps")

# Load lookup table based on species
if(species == "trichuris"){
  lookup_file <- file.path(kPathToMapsArtefacts, "iu_task_lookup_trichuris.rds")
}else if (species %in% c("ascaris","hookworm")){
  lookup_file <- file.path(kPathToMapsArtefacts, "iu_task_lookup_sth.rds")
}else if (species == "haematobium"){
  lookup_file <- file.path(kPathToMapsArtefacts, "iu_task_lookup_haema.rds")
} else if (species %in% c("mansoni_low_burden", "mansoni_high_burden")) {
  lookup_file <- file.path(kPathToMapsArtefacts, paste0("iu_task_lookup_", species, ".rds"))
} else {
  stop(paste("Unsupported species:", species))
}

if (!file.exists(lookup_file)) {
  stop(paste("Lookup file not found:", lookup_file))
}

load(lookup_file)
num_batches <- max(iu_task_lookup$TaskID)

cat(paste0("Total batches: ", num_batches, "\n"))

# Get list of batch IDs to analyze (excluding failed ones)
ids_to_analyze <- setdiff(1:num_batches, failed_ids)
cat(paste0("Analyzing ", length(ids_to_analyze), " batches\n\n"))

# Collect ESS data
ess_all_ius <- c()
missing_files <- c()
successful_batches <- c()

for(id in ids_to_analyze){
  
  # Construct filename based on our pipeline output format
  amis_file <- file.path(kPathToFittingArtefacts, paste0("fit_amis_", species, "_", id, "_sigma", sigma_value, ".RData"))
  
  if (!file.exists(amis_file)) {
    missing_files <- c(missing_files, id)
    next
  }
  
  tryCatch({
    load(amis_file) # loads amis_output
    
    if(exists("amis_output") && !is.null(amis_output$prevalence_map) && length(amis_output$prevalence_map) > 0) {
      iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
      ess <- amis_output$ess
      
      if(length(iu_names) == length(ess)) {
        ess_all_ius <- rbind(ess_all_ius, cbind(IU_ID = iu_names, TaskID = id, ESS = ess))
        successful_batches <- c(successful_batches, id)
      } else {
        cat(paste0("Warning: Mismatch in IU count vs ESS count for batch ", id, "\n"))
      }
    } else {
      cat(paste0("Warning: Invalid amis_output structure in batch ", id, "\n"))
    }
  }, error = function(e) {
    cat(paste0("Error loading batch ", id, ": ", e$message, "\n"))
  })
}

if(length(missing_files) > 0) {
  cat(paste0("Warning: AMIS output files not found for ", length(missing_files), " batches: ", paste(missing_files, collapse=", "), "\n"))
}

if(length(successful_batches) == 0) {
  stop("No valid AMIS output files found. Check that fitting has been completed and file paths are correct.")
}

cat(paste0("Successfully analyzed ", length(successful_batches), " batches\n\n"))

# Convert to data frame and analyze
ess_all_ius <- as.data.frame(ess_all_ius) %>%
  mutate_if(is.character, as.numeric)

# Find IUs with insufficient ESS
ess_ius_below_threshold <- ess_all_ius %>%
  filter(ESS < ess_threshold)

# Summary statistics
total_ius <- nrow(ess_all_ius)
low_ess_ius <- nrow(ess_ius_below_threshold)
low_ess_batches <- unique(ess_ius_below_threshold$TaskID)

cat("=== ESS ANALYSIS RESULTS ===\n")
cat(paste0("Total IUs analyzed: ", total_ius, "\n"))
cat(paste0("IUs with ESS < ", ess_threshold, ": ", low_ess_ius, " (", round(100*low_ess_ius/total_ius, 1), "%)\n"))
cat(paste0("Batches with insufficient ESS: ", length(low_ess_batches), "\n\n"))

if(length(low_ess_batches) > 0) {
  cat("BATCHES TO RERUN WITH HIGHER SIGMA:\n")
  cat(paste0(low_ess_batches, collapse=","), "\n\n")
  
  # Show detailed breakdown by batch
  cat("DETAILED BREAKDOWN:\n")
  batch_summary <- ess_ius_below_threshold %>%
    group_by(TaskID) %>%
    summarise(
      num_low_ess_ius = n(),
      min_ess = min(ESS),
      max_ess = max(ESS),
      avg_ess = round(mean(ESS), 1)
    ) %>%
    arrange(TaskID)
  
  for(i in 1:nrow(batch_summary)) {
    row <- batch_summary[i,]
    cat(paste0("Batch ", row$TaskID, ": ", row$num_low_ess_ius, " IUs with low ESS (min: ", 
               round(row$min_ess, 1), ", max: ", round(row$max_ess, 1), ", avg: ", row$avg_ess, ")\n"))
  }
} else {
  cat("✓ All batches have sufficient ESS!\n")
}

cat("\nRECOMMENDED ACTION:\n")
if(length(low_ess_batches) > 0) {
  cat("Rerun fitting for the flagged batches with higher sigma (e.g., 0.025):\n")
  cat("--stage=fitting --id=", paste(low_ess_batches, collapse=" --id="), " --amis-sigma=0.025\n")
} else {
  cat("No action needed - all batches meet the ESS threshold.\n")
}

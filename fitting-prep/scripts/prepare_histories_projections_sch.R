library(dplyr)
library(tidyr)
library(writexl)
library(optparse)
library(sf)

# Define command line options
option_list <- list(
  make_option(c("-s", "--species"),
    type = "character",
    default = "haematobium",
    help = "Species to process (haematobium, mansoni_low_burden, mansoni_high_burden) [default=%default]"
  ),
  make_option(c("-i", "--id"),
    type = "integer",
    default = NULL,
    help = "Single batch ID to process. If not provided, will process all batches"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Set up paths from environment variables
kPathToInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToMapsSCH <- file.path(kPathToInputs, "Maps-SCH")
kPathToMapsArtefacts <- file.path(kPathToArtefacts, "Maps")
kPathToESPEN <- file.path(kPathToInputs, "ESPEN_IU_2021")

# Create output directories if they don't exist
if (!dir.exists(kPathToMapsArtefacts)) {
  dir.create(kPathToMapsArtefacts, recursive = TRUE)
}

# Process species from command line argument - map to internal naming
species_input <- opts$species
if (species_input == "haematobium") {
  species_sch <- c("haematobium")
} else if (species_input %in% c("mansoni_low_burden", "mansoni_high_burden")) {
  species_sch <- c("mansoni")  # Both variants use same mansoni processing
} else {
  stop(paste("Invalid species:", species_input, ". Valid choices: haematobium, mansoni_low_burden, mansoni_high_burden"))
}

for (species in species_sch){
  if (species=="haematobium"){
    load(file.path(kPathToMapsArtefacts, "haematobium_maps.rds"))
    prevalence_map = get(paste0("haematobium_maps"))
  } else {
    load(file.path(kPathToMapsArtefacts, "mansoni_maps.rds"))
    prevalence_map = get(paste0("mansoni_maps"))
  }
  
  map_2002 = prevalence_map[[1]]$data
  
  # Read in MDA history
  sch_histories_raw = read.csv(file.path(kPathToMapsSCH, "Schisto_IU_Cleaned_1.csv")) %>%
    mutate(EpiCov_binned = case_when(
      EpiCov <= 15 ~ 0,
      EpiCov > 15 & EpiCov <= 75 ~ 0.15,
      EpiCov > 75 ~ 0.75
    ))
  sch_histories = sch_histories_raw %>%
    filter(IU_ID_MAPPING %in% map_2002$IU_ID)  %>%
    mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
    select(IU_ID_MAPPING,Year,TargetPop,Year_TargetPop,EpiCov_binned,MDA_scheme) 
  
  # backfill pre-2013 histories
  # data from: https://www.nejm.org/doi/suppl/10.1056/NEJMoa1812165/suppl_file/nejmoa1812165_appendix.pdf
  pre2013_countries_list = c("Burkina Faso","Burundi","Malawi","Mali-Segou","Mali-Bamako/Koulikoro","Niger","Rwanda","Tanzania","Uganda")
  pre2013_reps = data.frame(
    country = pre2013_countries_list,
    start_year = c(2004,2008,2012,2004,2004,2004,2008,2005,2003)) %>%
    mutate(reps = 2013-start_year)
  
  # take average coverage from 2013:2015 to backfill
  pre2013_ius = sch_histories_raw %>%
    filter(ADMIN0 %in% c("Burkina Faso","Burundi","Malawi","Mali","Niger","Rwanda","Tanzania (Mainland)","Tanzania (Zanzibar)","Uganda")) %>%
    filter(IU_ID_MAPPING %in% map_2002$IU_ID)  %>%
    filter(Year %in% 2013:2015 & PC == 1) %>%
    filter((!ADMIN0 %in% c("Mali")) | (ADMIN0 == "Mali" & ADMIN1 %in% c("Segou","Bamako","Koulikoro"))) %>%
    mutate(country = ifelse(ADMIN0 %in% c("Tanzania (Mainland)","Tanzania (Zanzibar)"), "Tanzania",
                            ifelse(ADMIN0 == "Mali" & ADMIN1 == "Segou", "Mali-Segou",
                                   ifelse(ADMIN0 == "Mali" & ADMIN1 %in% c("Bamako","Koulikoro"), "Mali-Bamako/Koulikoro", ADMIN0)))) %>%
    group_by(country, ADMIN1, IU_ID_MAPPING, TargetPop, MDA_scheme) %>%
    summarise(mean_EpiCov = mean(EpiCov)) %>%
    ungroup() %>%
    mutate(EpiCov_binned = case_when(
      mean_EpiCov <= 15 ~ 0,
      mean_EpiCov > 15 & mean_EpiCov <= 75 ~ 0.15,
      mean_EpiCov > 75 ~ 0.75
    )) %>%
    select(country, ADMIN1, IU_ID_MAPPING, TargetPop, EpiCov_binned, MDA_scheme)
  
  # duplicate rows for Burundi Pilot
  # 3 provinces: Ortu G, Assoum M, Wittmann U, Knowles S, Clements M, Ndayishimiye O, Basáñez MG, Lau C, Clements A, Fenwick A, Magalhaes RJ. The impact of an 8-year mass drug administration programme on prevalence, intensity and co-infections of soil-transmitted helminthiases in Burundi. Parasit Vectors. 2016 Sep 22;9(1):513. doi: 10.1186/s13071-016-1794-9. PMID: 27660114; PMCID: PMC5034474.
  # (not used) coverage = 52% in Ndayishimiye O, Ortu G, Soares Magalhaes RJ, Clements A, Willems J, et al. (2014) Control of Neglected Tropical Diseases in Burundi: Partnerships, Achievements, Challenges, and Lessons Learned after Four Years of Programme Implementation. PLOS Neglected Tropical Diseases 8(5): e2684. https://doi.org/10.1371/journal.pntd.0002684
  burundi_pilot_2007 = pre2013_ius %>%
    filter(IU_ID_MAPPING %in% unique(pre2013_ius %>% 
                                       filter(ADMIN1 %in% c("Bururi","Bubanza","Cibitoke")) %>%
                                       select(IU_ID_MAPPING))$IU_ID_MAPPING) %>%
    mutate(Year=2007) %>%
    select(country, Year, ADMIN1, IU_ID_MAPPING, TargetPop, EpiCov_binned, MDA_scheme)
  
  # join datasets
  # note: no map samples for Niger and Uganda
  pre2013_mda_data = rbind(
    data.frame(country = rep(pre2013_countries_list, pre2013_reps$reps),
               Year = c(2004:2012,2008:2012,2012,2004:2012,2004:2012,2004:2012,2008:2012,2005:2012,2003:2012)) %>%
      left_join(pre2013_ius, by=c("country"), relationship = "many-to-many"), 
    burundi_pilot_2007) %>%
    filter((!is.na(IU_ID_MAPPING)) & TargetPop %in% c("SAC","Adults")) %>%
    mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
    select(IU_ID_MAPPING, Year, TargetPop, Year_TargetPop, EpiCov_binned, MDA_scheme)
  
  # Combine both histories
  histories_joined = rbind(sch_histories,pre2013_mda_data) %>%
    arrange(IU_ID_MAPPING,Year)
  histories_joined$MDA_scheme[histories_joined$EpiCov_binned ==0 & !histories_joined$MDA_scheme=="Not delivered"] = "Not delivered"
  histories_joined$EpiCov_binned[histories_joined$MDA_scheme=="Not delivered"] = NA
  
  # group drugs 
  histories_joined$MDA_scheme[histories_joined$MDA_scheme == "PZQ+ALB/MBD"] = "PZQ"
  histories_joined$MDA_scheme[!histories_joined$MDA_scheme %in% c("PZQ","Not delivered")] = "notPZQ"
  # rename drugs
  PZQ_name = "Old Product B (SOC)" 
  notPZQ_name = "New Product A" 
  histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("PZQ")] = PZQ_name
  histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("notPZQ")] = notPZQ_name
  
  add_burnin_year = data.frame(IU_ID_MAPPING = map_2002$IU_ID,
                               Year = 1985,
                               TargetPop = "Adults",
                               Year_TargetPop = "1985_Adults",
                               EpiCov_binned = NA,
                               MDA_scheme = "Start burn in")
  histories_joined = rbind(histories_joined,add_burnin_year) %>%
    arrange(IU_ID_MAPPING,Year)
  
  # Repeat most recent history if any treatment was given in the previous 3 years
  recent_treatments_since_2020 = histories_joined %>%
    filter(!MDA_scheme == "Not delivered") %>%
    group_by(IU_ID_MAPPING) %>%
    summarise(last_treated_year = floor(max(Year))) %>%
    filter(last_treated_year >= 2020) %>%
    left_join(histories_joined %>% mutate(floor_year = floor(Year)), by=c("IU_ID_MAPPING","last_treated_year"="floor_year")) %>%
    select(IU_ID_MAPPING, TargetPop, Year_TargetPop, EpiCov_binned, MDA_scheme)
  recent_treatments_sch = recent_treatments_since_2020 %>% filter(TargetPop %in% c("SAC","Adults"))
  
  recent_treatments_future_years = cbind(Year = c(rep(2023:2025, each=nrow(recent_treatments_sch))),
                                         rbind(recent_treatments_sch,recent_treatments_sch,recent_treatments_sch)) %>%
    mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
    select(IU_ID_MAPPING, Year, TargetPop, Year_TargetPop, EpiCov_binned, MDA_scheme)
  
  # Join data
  projections_data = rbind(histories_joined, recent_treatments_future_years)
  projections_full = data.frame(IU_ID_MAPPING = rep(map_2002$IU_ID,each=47),
                                Year = rep(c(1985,rep(2003:2025,each=2)), length(map_2002$IU_ID)),
                                TargetPop = rep(c("Adults",rep(c("SAC", "Adults"),23)),length(map_2002$IU_ID))) %>%
    mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>% 
    left_join(projections_data) 
  projections_full$MDA_scheme[is.na(projections_full$MDA_scheme)] = "Not delivered"
  
  # Load task IDs and Overwrite MDA files for IUs with no treatments to add negligible coverage (otherwise python code will break)
  if (species == "haematobium"){
    load(file.path(kPathToMapsArtefacts, "iu_task_lookup_haema.rds"))
    id_no_mda = c(1124,1227)
  } else {
    # For mansoni, use the species_input which has the variant (mansoni_low_burden or mansoni_high_burden)
    lookup_file <- paste0("iu_task_lookup_", species_input, ".rds")
    load(file.path(kPathToMapsArtefacts, lookup_file))
    id_no_mda = c(1088,1204)
  }
  
  iu_no_mda = unique((projections_full %>% 
                        filter(IU_ID_MAPPING %in% (iu_task_lookup %>% 
                                                     filter(TaskID %in% id_no_mda))$IU_ID))$IU_ID_MAPPING)
  inds = which(projections_full$IU_ID_MAPPING %in% iu_no_mda & projections_full$Year %in% c(2003,2004))
  projections_full[inds, colnames(projections_full) %in% c("EpiCov_binned")] = 1e-11
  projections_full[inds, colnames(projections_full) %in% c("MDA_scheme")] = notPZQ_name
  
  write.csv(projections_full, file=file.path(kPathToMapsArtefacts, paste0("mda_history_",species,".csv")))
  
  # Get into format for python model
  coverage_wide = projections_full %>%
    filter(!MDA_scheme == "Not delivered") %>%
    select(-c(Year_TargetPop,MDA_scheme)) %>%
    arrange(Year) %>%
    pivot_wider(names_from = Year , values_from = EpiCov_binned) %>%
    mutate("Country/Region" = "All",
           "Intervention Type"= "Treatment",
           "Platform Type" = "Campaign",
           Platform = "MDA",
           #Drug = MDA_scheme,
           "Cohort (if not total pop in country/region)" = NA,
           "min age" = case_when(
             TargetPop == "Adults" ~ 15,
             TargetPop == "SAC" ~ 5
           ),
           "max age" = case_when(
             TargetPop == "Adults" ~ 100,
             TargetPop == "SAC" ~ 15
           )) %>%
    select(IU_ID_MAPPING,"Country/Region","Intervention Type","Platform Type", Platform, #Drug,
           "Cohort (if not total pop in country/region)","min age","max age",
           starts_with(c("19","20")))
  
  # Get into format for python model
  drug_wide = unique(projections_full %>%
                       filter(!MDA_scheme == "Not delivered") %>%
                       select(IU_ID_MAPPING,Year,MDA_scheme)) %>%
    arrange(Year) %>%
    mutate(drug_indicator = ifelse(MDA_scheme == "Not delivered", NA,1)) %>%
    pivot_wider(names_from = Year , values_from = drug_indicator) %>%
    mutate("Country/Region" = "All",
           Platform = "MDA",
           Product = MDA_scheme) %>%
    select(IU_ID_MAPPING,"Country/Region", Platform, Product, starts_with(c("19","20"))) %>%
    filter(!Product == "Start burn in")
  
  
  # export histories files and define batches to run
  if (species == "haematobium"){
    inputs_path = file.path(kPathToArtefacts, "endgame_inputs", "sch-haematobium")
  } else {
    inputs_path = file.path(kPathToArtefacts, "endgame_inputs", "sch-mansoni")
  }
  if (!dir.exists(inputs_path)) {dir.create(inputs_path, recursive = TRUE)}
  
  # Determine which batches to process
  if (!is.null(opts$id)) {
    # Validate single batch ID
    max_batch_id <- max(iu_task_lookup$TaskID)
    if (opts$id > max_batch_id || opts$id < 1) {
      stop(paste("Specified batch ID", opts$id, "is out of range. Valid range: 1 to", max_batch_id))
    }
    id_list <- opts$id
    cat(paste0("Processing single batch ID: ", opts$id, " for species: ", species, "\n"))
  } else {
    id_list <- 1:max(iu_task_lookup$TaskID)
    cat(paste0("Processing all ", max(iu_task_lookup$TaskID), " batches for species: ", species, "\n"))
  }
  
  
  for (id in id_list){
    
    ius_per_batch = iu_task_lookup %>%
      filter(TaskID == id)
    
    ius = matrix(ius_per_batch$IU_ID,ncol=1)
    
    for (iu in ius) {
      
      mda_cov_iu = coverage_wide %>%
        filter(IU_ID_MAPPING == iu) %>%
        select(-IU_ID_MAPPING) 
      
      mda_cov_iu_xl = rbind(colnames(mda_cov_iu),mda_cov_iu) %>% 
        mutate_at(vars(starts_with(c("19","20"))), as.numeric)
      
      drug_template = data.frame("Country/Region"	= "All",
                                 Platform	= "MDA",
                                 Product = c("Old Product B (SOC)","New Product A"), check.names = F)
      
      mda_drug_iu = suppressMessages( drug_template %>%
                                        left_join((drug_wide %>%
                                                     filter(IU_ID_MAPPING == iu) %>%
                                                     select(-IU_ID_MAPPING) )))
      
      mda_drug_iu_xl = rbind(colnames(mda_drug_iu),mda_drug_iu) %>% 
        mutate_at(vars(starts_with(c("19","20"))), as.numeric)
      
      # The first year columns must be a numeric integer for it to be read correctly into python!!
      mda_path<-file.path(inputs_path,paste0("InputMDA_MTP_projections_",iu,".xlsx"))
      write_xlsx(list("Platform Coverage" = mda_cov_iu_xl, "MarketShare" = mda_drug_iu_xl),path=mda_path, col_names = F) # write input MDA file
    }
  }
  
  # make look up table for country codes
  table_iu_idx <- prevalence_map[[1]]$data[,c("IU_ID","TaskID")]
  colnames(table_iu_idx) <- c("IU_CODE","TaskID")
  rownames(table_iu_idx) <- NULL
  df <- read_sf(dsn = kPathToESPEN, layer = "ESPEN_IU_2021") %>%
    filter(IU_ID %in% table_iu_idx$IU_CODE)
  df <- st_drop_geometry(df[,c("IU_ID","ADMIN0ISO3")])
  head(df)
  table_iu_idx$country <- NA
  for(i in 1:nrow(table_iu_idx)){
    IU <- table_iu_idx[i,"IU_CODE"]
    wh <- which(df$IU_ID==IU)[1]
    country <- df[wh,"ADMIN0ISO3"]
    table_iu_idx[i,"country"] <- country
  }
  table_iu_idx$IU_CODE <- as.character(table_iu_idx$IU_CODE)
  table_iu_idx$TaskID <- as.integer(table_iu_idx$TaskID)
  table_iu_idx$country <- as.character(table_iu_idx$country)
  
  write.csv(table_iu_idx, file=file.path(kPathToMapsArtefacts, paste0("table_iu_idx_", species, ".csv")), row.names=F)
  
}

cat("Finished running prepare_histories_projections_sch.R \n")

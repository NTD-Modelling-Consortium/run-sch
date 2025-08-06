library(dplyr)
library(tidyr)
library(writexl)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-s", "--species"),
    type = "character",
    default = "ascaris",
    help = "Species to process (ascaris, hookworm) [default=%default]"
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
kPathToMapsSTH <- file.path(kPathToInputs, "Maps-STH")
kPathToMapsArtefacts <- file.path(kPathToArtefacts, "Maps")

# Create output directory if it doesn't exist
if (!dir.exists(kPathToMapsArtefacts)) {
  dir.create(kPathToMapsArtefacts, recursive = TRUE)
}

# Read in data and remove problematic IUs
sth_ius = unique(read.csv(file.path(kPathToMapsSTH, "sartorious_2021.csv")) %>%
  filter(!p4==1) %>%
  select(IU_2021))

sth_histories_raw = read.csv(file.path(kPathToMapsSTH, "STHCleaned_1.csv")) 
sth_histories = sth_histories_raw %>%
  filter(IU_ID_MAPPING %in% sth_ius$IU_2021) %>%
  mutate(EpiCov_binned = case_when(
    EpiCov <= 15 ~ 0,
    EpiCov > 15 & EpiCov <= 75 ~ 0.15,
    EpiCov > 75 ~ 0.75
  )) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
  select(IU_ID_MAPPING,Year,TargetPop,Year_TargetPop,EpiCov_binned,MDA_scheme) 

# Use LF histories for pre2013
lf_histories_raw= read.csv(file.path(kPathToMapsSTH, "LF_MDA_Africa_2024_IU_updated.csv")) 
lf_histories = lf_histories_raw %>%
  mutate(IU_ID_MAPPING = as.numeric(substr(IU_ID,4,8))) %>%
  mutate_at(vars(starts_with("Cov")), ~cut(., 
                                          breaks=c(-Inf, 15, 75, Inf), 
                                          labels=c(0 ,0.15, 0.75))) %>%
  select(-c(IU_ID,ADMIN0,IUs_NAME)) %>%
  filter(IU_ID_MAPPING %in% sth_ius$IU_2021)
  
lf_cov_long = lf_histories %>%
  select(IU_ID_MAPPING,starts_with("Cov")) %>%
  pivot_longer(starts_with("Cov"), values_to = "EpiCov_binned") %>%
  mutate(Year = as.numeric(substr(name,4,7)) + 0.5) %>%
  select(IU_ID_MAPPING, Year, EpiCov_binned)

lf_MDA_long = lf_histories %>%
  select(IU_ID_MAPPING,starts_with("MDA")) %>%
  mutate(MDA_scheme2000 = MDA_schemebefore2013,
         MDA_scheme2001 = MDA_schemebefore2013,
         MDA_scheme2002 = MDA_schemebefore2013,
         MDA_scheme2003 = MDA_schemebefore2013,
         MDA_scheme2004 = MDA_schemebefore2013,
         MDA_scheme2005 = MDA_schemebefore2013,
         MDA_scheme2006 = MDA_schemebefore2013,
         MDA_scheme2007 = MDA_schemebefore2013,
         MDA_scheme2008 = MDA_schemebefore2013,
         MDA_scheme2009 = MDA_schemebefore2013,
         MDA_scheme2010 = MDA_schemebefore2013,
         MDA_scheme2011 = MDA_schemebefore2013,
         MDA_scheme2012 = MDA_schemebefore2013) %>%
  select(-MDA_schemebefore2013) %>%
  pivot_longer(starts_with("MDA"), values_to = "MDA_scheme") %>%
  mutate(Year = as.numeric(substr(name,11,nchar(name))) + 0.5) %>%
  select(IU_ID_MAPPING, Year, MDA_scheme)

lf_joined = lf_cov_long %>%
  mutate(TargetPop = "SAC/Adults") %>% 
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
  select(IU_ID_MAPPING, Year, TargetPop, Year_TargetPop, EpiCov_binned) %>%
  left_join(lf_MDA_long) 

# Combine both histories
histories_joined =  rbind(sth_histories, lf_joined) %>%
  arrange(IU_ID_MAPPING,Year)
histories_joined$MDA_scheme[histories_joined$EpiCov_binned ==0 & !histories_joined$MDA_scheme=="Not delivered"] = "Not delivered"
histories_joined$EpiCov_binned[histories_joined$MDA_scheme=="Not delivered"] = NA
histories_joined$MDA_scheme[histories_joined$MDA_scheme == "PZQ+ALB/MBD"] = "ALB/MBD"
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("ALB+DEC")] = "ALB"

# select based on disease for now
# for ascaris and hookworm
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("IDA","ALB+IVM")] = "ALB" 
# rename drugs
ALB_name = "Old Product B (SOC)" 
MBD_name = "New Product A" 
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("ALB")] = ALB_name
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("ALB/MBD")] = MBD_name

add_burnin_year = data.frame(IU_ID_MAPPING = sth_ius$IU_2021,
                        Year = 1985,
                        TargetPop = "SAC/Adults",
                        Year_TargetPop = "1985_SAC/Adults",
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
recent_treatments_sch = recent_treatments_since_2020 %>% filter(TargetPop %in% c("PreSAC","SAC"))
recent_treatments_lf = recent_treatments_since_2020 %>% filter(TargetPop %in% c("SAC/Adults"))

recent_treatments_future_years = cbind(Year = c(rep(2023:2025, each=nrow(recent_treatments_sch)), 
                                                rep(seq(2023.5,2025.5,by=1), each=nrow(recent_treatments_lf))),
                                       rbind(recent_treatments_sch,recent_treatments_sch,recent_treatments_sch,
                                             recent_treatments_lf,recent_treatments_lf,recent_treatments_lf)) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
  select(IU_ID_MAPPING, Year, TargetPop, Year_TargetPop, EpiCov_binned, MDA_scheme)

# Join data
projections_data = rbind(histories_joined, recent_treatments_future_years)
projections_full = data.frame(IU_ID_MAPPING = rep(sth_ius$IU_2021,each=66),
                              Year = rep(c(1985,seq(2000,2012.5,by=0.5),rep(2013:2025,each=2),seq(2013.5,2025.5,by=1)), length(sth_ius$IU_2021)),
                              TargetPop = c(rep("SAC/Adults",27),rep(c("PreSAC", "SAC"),13),rep("SAC/Adults",13))) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>% 
  left_join(projections_data) 
projections_full$MDA_scheme[is.na(projections_full$MDA_scheme)] = "Not delivered"

# Load task IDs
load(file.path(kPathToMapsArtefacts, "iu_task_lookup_sth.rds"))

# Overwrite MDA files for IUs with no treatments to add negligible coverage (otherwise python code will break)
iu_no_mda = unique((projections_full %>% 
                      filter(IU_ID_MAPPING %in% (iu_task_lookup %>% 
                                                   filter(TaskID == max(iu_task_lookup$TaskID)))$IU_2021))$IU_ID_MAPPING)
inds = which(projections_full$IU_ID_MAPPING %in% iu_no_mda & projections_full$Year %in% c(2000.5,2001.5))
projections_full[inds, colnames(projections_full) %in% c("EpiCov_binned")] = 1e-11
projections_full[inds, colnames(projections_full) %in% c("MDA_scheme")] = MBD_name

# Duplicate row when MDA_scheme = ALBx2
albx2_rows = projections_full %>%
  filter(MDA_scheme == "ALBx2") %>%
  mutate(Year = Year - 0.49) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
  select(IU_ID_MAPPING,Year,TargetPop,Year_TargetPop,EpiCov_binned,MDA_scheme)
projections_full = rbind(projections_full,albx2_rows) %>%
  mutate(MDA_scheme = ifelse(MDA_scheme=="ALBx2",ALB_name,MDA_scheme))

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
           TargetPop == "SAC/Adults" ~ 5,
           TargetPop == "PreSAC" ~ 0,
           TargetPop == "SAC" ~ 5
         ),
         "max age" = case_when(
           TargetPop == "SAC/Adults" ~ 100,
           TargetPop == "PreSAC" ~ 5,
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


# export histories files
inputs_path = file.path(kPathToArtefacts, "endgame_inputs", "STH")
if (!dir.exists(inputs_path)) {dir.create(inputs_path, recursive = TRUE)}

# Define new batches for projections   # iu_task_lookup was the original batch-IU lookup table used for fitting
IUs <- sort(unique(iu_task_lookup$IU_2021))
num_IUs <- length(IUs)
batch_size <- 1 # number of IUs in each batch for projections
num_batches <- num_IUs/batch_size
proj_iu_task_lookup <- data.frame(IU_2021=IUs, TaskID=sort(rep(1:num_batches, batch_size)))
save(proj_iu_task_lookup, file=file.path(kPathToMapsArtefacts, "proj_iu_task_lookup_STH.rds"))

num_batches <- max(proj_iu_task_lookup$TaskID)
cat(paste0("Number of batches for projections: ", num_batches, "\n"))

# Determine which batches to process
if (!is.null(opts$id)) {
  # Validate single batch ID
  if (opts$id > num_batches || opts$id < 1) {
    stop(paste("Specified batch ID", opts$id, "is out of range. Valid range: 1 to", num_batches))
  }
  batch_ids <- opts$id
  cat(paste0("Processing single batch ID: ", opts$id, "\n"))
} else {
  batch_ids <- 1:num_batches
  cat(paste0("Processing all ", num_batches, " batches\n"))
}

for (id in batch_ids){

  ius_per_batch = proj_iu_task_lookup %>%
    filter(TaskID == id)
  ius = matrix(ius_per_batch$IU_2021,ncol=1)
  
  for (iu in ius[,1]) {
    
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
  cat(paste0("Files prepared for id = ", id, "; "))
}

# Use species from command line argument
species <- opts$species

# Validate species choice
valid_species <- c("ascaris", "hookworm", "trichuris")
if (!species %in% valid_species) {
  stop(paste("Invalid species:", species, ". Valid choices:", paste(valid_species, collapse=", ")))
}

# read in map
load(file.path(kPathToMapsArtefacts, paste0(species,'_maps.rds'))) # load species_map_allyears

if(species=="ascaris"){
  species_map_allyears <- ascaris_map_allyears
}else if(species=="hookworm"){
  species_map_allyears <- hookworm_map_allyears
}else if(species=="trichuris"){
  species_map_allyears <- trichuris_map_allyears
}

# Look up table that links each IU to its index in the corresponding (fitting) batch
table_iu_idx <- species_map_allyears[[1]]$data[,c("IU_2021","TaskID")]
colnames(table_iu_idx) <- c("IU_CODE","TaskID")
rownames(table_iu_idx) <- NULL
df <- read.csv(file.path(kPathToMapsSTH, "sartorious_2021.csv"))
df <- df[,c("IU_2021","ADMIN0ISO3")]
head(df)
table_iu_idx$country <- NA
for(i in 1:nrow(table_iu_idx)){
  IU <- table_iu_idx[i,"IU_CODE"]
  wh <- which(df$IU_2021==IU)[1]
  country <- df[wh,"ADMIN0ISO3"]
  table_iu_idx[i,"country"] <- country
}
table_iu_idx$IU_CODE <- as.character(table_iu_idx$IU_CODE)
table_iu_idx$TaskID <- as.integer(table_iu_idx$TaskID)
table_iu_idx$country <- as.character(table_iu_idx$country)
write.csv(table_iu_idx, file = file.path(kPathToMapsArtefacts, "table_iu_idx_STH.csv"), row.names = F)

cat("Finished running prepare_histories_projections.R \n")


library(dplyr)
library(tidyr)
library(writexl)

# Get paths from environment variables
kPathToInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToMapsSTH <- file.path(kPathToInputs, "Maps-STH")
kPathToMapsArtefacts <- file.path(kPathToArtefacts, "Maps")

# Create output directories if they don't exist
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
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("PZQ+ALB/MBD","ALB+DEC","ALB/MBD")] = "ALB"
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("IDA")] = "ALB+IVM"
# rename drugs
ALB_name = "Old Product B (SOC)" 
IVM_name = "New Product A" 
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("ALB")] = ALB_name
histories_joined$MDA_scheme[histories_joined$MDA_scheme %in% c("ALB+IVM")] = IVM_name

add_burnin_year = data.frame(IU_ID_MAPPING = sth_ius$IU_2021,
                        Year = 1985,
                        TargetPop = "SAC/Adults",
                        Year_TargetPop = "1985_SAC/Adults",
                        EpiCov_binned = NA,
                        MDA_scheme = "Start burn in")
histories_joined = rbind(histories_joined,add_burnin_year) %>%
  arrange(IU_ID_MAPPING,Year)


histories_joined_full = data.frame(IU_ID_MAPPING = rep(sth_ius$IU_2021,each=57),
                              Year = rep(c(1985,seq(2000,2012.5,by=0.5),rep(2013:2022,each=2),seq(2013.5,2022.5,by=1)), length(sth_ius$IU_2021)),
                              TargetPop = c(rep("SAC/Adults",27),rep(c("PreSAC", "SAC"),10),rep("SAC/Adults",10))) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>% 
  left_join(histories_joined) 
histories_joined_full$MDA_scheme[is.na(histories_joined_full$MDA_scheme)] = "Not delivered"


# Determine IUs with the same treatment history to assign batches

# find unique coverage histories
coverage_history = histories_joined_full %>%
  filter(Year < 2019) %>%
  select(IU_ID_MAPPING,Year_TargetPop,EpiCov_binned) %>%
  pivot_wider(names_from="Year_TargetPop",values_from="EpiCov_binned")

unique_coverage = coverage_history %>% 
  group_by(across(c(-IU_ID_MAPPING))) %>% 
  summarise(ius= list(IU_ID_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

coverage_indexes = lapply(1:nrow(unique_coverage), function(i) data.frame(unique_coverage[["ius"]][[i]],unique_coverage[["index"]][i]))
coverage_indexes = do.call(rbind,coverage_indexes)
colnames(coverage_indexes) = c("IU","coverage_index")


# find unique MDA_scheme histories
MDA_scheme_history = histories_joined_full %>%
  filter(Year < 2019) %>%
  select(IU_ID_MAPPING,Year_TargetPop,MDA_scheme) %>%
  pivot_wider(names_from="Year_TargetPop",values_from="MDA_scheme")

unique_MDA_scheme = MDA_scheme_history %>% 
  group_by(across(c(-IU_ID_MAPPING))) %>% 
  summarise(ius= list(IU_ID_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

MDA_scheme_indexes = lapply(1:nrow(unique_MDA_scheme), function(i) data.frame(unique_MDA_scheme[["ius"]][[i]],unique_MDA_scheme[["index"]][i]))
MDA_scheme_indexes = do.call(rbind,MDA_scheme_indexes)
colnames(MDA_scheme_indexes) = c("IU","MDA_scheme_index")


# join all and find unique combinations
all_treatments = data.frame(IU=sth_ius$IU_2021) %>%
  left_join(coverage_indexes) %>%
  left_join(MDA_scheme_indexes)

unique_all = all_treatments %>%
  group_by(across(c(-IU))) %>% 
  summarise(ius= list(IU)) %>%
  ungroup() %>%
  mutate(index = row_number())

# lookup table for IUs to batch IDs
iu_task_lookup = lapply(1:nrow(unique_all), function(i) data.frame(unique_all[["ius"]][[i]],unique_all[["index"]][i]))
iu_task_lookup = do.call(rbind,iu_task_lookup) 
colnames(iu_task_lookup) = c("IU_2021","TaskID")
# re-define batches for IUs with very high prevalances and 
# and tight intervals in the year 2000 ---------------------
original_last_batch_ID <- max(iu_task_lookup$TaskID)
wh <- which(iu_task_lookup$IU_2021==19448)
iu_task_lookup[wh,"TaskID"] <- original_last_batch_ID + 1L
wh <- which(iu_task_lookup$IU_2021==19517)
iu_task_lookup[wh,"TaskID"] <- original_last_batch_ID + 2L
# ----------------------------------------------------------
iu_task_lookup = iu_task_lookup %>%
  left_join(unique(sth_histories %>% select(IU_ID_MAPPING)), by=c("IU_2021"="IU_ID_MAPPING"))
save(iu_task_lookup, file=file.path(kPathToMapsArtefacts, "iu_task_lookup_trichuris.rds"))

# save batch IDs to maps file
sth_distributions = iu_task_lookup %>% left_join(read.csv(file.path(kPathToMapsSTH, "sartorious_2021.csv")) %>% filter(!p4==1))
write.csv(sth_distributions,file=file.path(kPathToMapsArtefacts, "sartorious_2021_withTaskID_trichuris.csv"), row.names = F)


# Overwrite MDA files for IUs with no treatments to add negligible coverage (otherwise python code will break)
iu_no_mda = unique((histories_joined_full %>% 
                      filter(IU_ID_MAPPING %in% (iu_task_lookup %>% 
                                                   filter(TaskID == original_last_batch_ID))$IU_2021))$IU_ID_MAPPING)
inds = which(histories_joined_full$IU_ID_MAPPING %in% iu_no_mda & histories_joined_full$Year %in% c(2000.5,2001.5))
histories_joined_full[inds, colnames(histories_joined_full) %in% c("EpiCov_binned")] = 1e-11
histories_joined_full[inds, colnames(histories_joined_full) %in% c("MDA_scheme")] = IVM_name

# Duplicate row when MDA_scheme = ALBx2
albx2_rows = histories_joined_full %>%
  filter(MDA_scheme == "ALBx2") %>%
  mutate(Year = Year - 0.49) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>%
  select(IU_ID_MAPPING,Year,TargetPop,Year_TargetPop,EpiCov_binned,MDA_scheme)
histories_joined_full = rbind(histories_joined_full,albx2_rows) %>%
  mutate(MDA_scheme = ifelse(MDA_scheme=="ALBx2",ALB_name,MDA_scheme))


# Get into format for python model
coverage_wide = histories_joined_full %>%
  filter(!MDA_scheme == "Not delivered") %>%
  select(-c(Year_TargetPop,MDA_scheme)) %>%
  filter(Year < 2019) %>%
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
drug_wide = unique(histories_joined_full %>%
  filter(!MDA_scheme == "Not delivered") %>%
  select(IU_ID_MAPPING,Year,MDA_scheme)) %>%
  filter(Year < 2019) %>%
  arrange(Year) %>%
  mutate(drug_indicator = ifelse(MDA_scheme == "Not delivered", NA,1)) %>%
  pivot_wider(names_from = Year , values_from = drug_indicator) %>%
  mutate("Country/Region" = "All",
         Platform = "MDA",
         Product = MDA_scheme) %>%
  select(IU_ID_MAPPING,"Country/Region", Platform, Product, starts_with(c("19","20"))) %>%
  filter(!Product == "Start burn in")


# export histories files
kPathToEndgameInputs <- file.path(kPathToArtefacts, "endgame_inputs", "STH")
if (!dir.exists(kPathToEndgameInputs)) {dir.create(kPathToEndgameInputs, recursive = TRUE)}

for (id in 1:max(iu_task_lookup$TaskID)){

  iu_file<-file.path(kPathToEndgameInputs,paste0("IUs_MTP_trichuris_",id,".csv"))
  ius_per_batch = iu_task_lookup %>%
    filter(TaskID == id)
  ius = matrix(ius_per_batch$IU_2021,ncol=1)
  write.table(ius,file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
  
  for (iu in ius[1,1]) { # all IUs in the same batch have the same treatment
    
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
    mda_path<-file.path(kPathToEndgameInputs,paste0("InputMDA_MTP_trichuris_",id,".xlsx"))
    write_xlsx(list("Platform Coverage" = mda_cov_iu_xl, "MarketShare" = mda_drug_iu_xl),path=mda_path, col_names = F) # write input MDA file
  }
}

cat("Finished running prepare_histories_trichuris.R \n")

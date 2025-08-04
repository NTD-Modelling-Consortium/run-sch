library(dplyr)
library(tidyr)
library(readr)
library(writexl)
library(readxl)
library(optparse)

# Command line arguments
option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    default = NULL,
    help = "Optional batch ID to process. If not provided, all batches will be processed."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Get paths from environment variables
kPathToInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToMapsSCH <- file.path(kPathToInputs, "Maps-SCH")
kPathToMapsArtefacts <- file.path(kPathToArtefacts, "Maps")

# Create output directories if they don't exist
if (!dir.exists(kPathToMapsArtefacts)) {
  dir.create(kPathToMapsArtefacts, recursive = TRUE)
}

############## Haematobium maps ############## 

# Combine 2023 maps
haematobium_g1_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G1_2023.xlsx"), col_types = "text") 
haematobium_g2_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G2_2023.xlsx"), col_types = "text") 
haematobium_g3_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G3_2023.xlsx"), col_types = "text") 
haematobium_g4_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G4_2023.xlsx"), col_types = "text") 
haematobium_g5_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G5_2023.xlsx"), col_types = "text") 
haematobium_g6_2023 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G6_2023.xlsx"), col_types = "text") 
colnames(haematobium_g1_2023) = colnames(haematobium_g2_2023)

# rename 2023 to 2022 because only MDA history up to 2022 was used
haematobium_2022 = rbind(haematobium_g1_2023,
                         haematobium_g2_2023,
                         haematobium_g3_2023,
                         haematobium_g4_2023,
                         haematobium_g5_2023,
                         haematobium_g6_2023)

haematobium_2022 = haematobium_2022[,3:ncol(haematobium_2022)]
colnames(haematobium_2022) = c("IU_ID",sprintf("X%s",1:1000))
haematobium_2022 =  haematobium_2022 %>%
  filter(!X1 %in% c("NA","NaN")) %>% 
  mutate_all(as.numeric) %>%
  mutate_at(vars(starts_with("X")), ~ ifelse(. > 0.99, NA, .))

# Combine 2005 maps
haematobium_g1_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G1_2005.xlsx"), col_types = "text") 
haematobium_g2_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G2_2005.xlsx"), col_types = "text") 
haematobium_g3_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G3_2005.xlsx"), col_types = "text") 
haematobium_g4_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G4_2005.xlsx"), col_types = "text") 
haematobium_g5_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G5_2005.xlsx"), col_types = "text") 
haematobium_g6_2005 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G6_2005.xlsx"), col_types = "text") 
colnames(haematobium_g1_2005) = colnames(haematobium_g2_2005)

# rename 2005 to 2003 as this was when the first MDA intervention occurred
haematobium_2002 = rbind(haematobium_g1_2005,
                         haematobium_g2_2005,
                         haematobium_g3_2005,
                         haematobium_g4_2005,
                         haematobium_g5_2005,
                         haematobium_g6_2005)

haematobium_2002 = haematobium_2002[,3:ncol(haematobium_2002)]
colnames(haematobium_2002) = c("IU_ID",sprintf("X%s",1:1000))
haematobium_2002 =  haematobium_2002 %>%
  filter(!X1 %in% c("NA","NaN")) %>% 
  mutate_all(as.numeric) 
haematobium_2002 = left_join(haematobium_2022 %>% select(IU_ID),haematobium_2002) %>%
  mutate_at(vars(starts_with("X")), ~ ifelse(. > 0.99, NA, .))


# Combine 2013 maps
haematobium_g1_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G1_2013.xlsx"), col_types = "text") 
haematobium_g2_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G2_2013.xlsx"), col_types = "text") 
haematobium_g3_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G3_2013.xlsx"), col_types = "text") 
haematobium_g4_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G4_2013.xlsx"), col_types = "text") 
haematobium_g5_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G5_2013.xlsx"), col_types = "text") 
haematobium_g6_2013 = read_excel(file.path(kPathToMapsSCH, "S.haematobium_G6_2013.xlsx"), col_types = "text") 
colnames(haematobium_g1_2013) = colnames(haematobium_g2_2013)

haematobium_2013 = rbind(haematobium_g1_2013,
                         haematobium_g2_2013,
                         haematobium_g3_2013,
                         haematobium_g4_2013,
                         haematobium_g5_2013,
                         haematobium_g6_2013)

haematobium_2013 = haematobium_2013[,3:ncol(haematobium_2013)]
colnames(haematobium_2013) = c("IU_ID",sprintf("X%s",1:1000))
haematobium_2013 =  haematobium_2013 %>%
  filter(!X1 %in% c("NA","NaN")) %>% 
  mutate_all(as.numeric) 
haematobium_2013 = left_join(haematobium_2022 %>% select(IU_ID),haematobium_2013) %>%
  mutate_at(vars(starts_with("X")), ~ ifelse(. > 0.99, NA, .))


##### identify IUs with only high prevalences (>0.99) that we can't fit
# summarise IU prevalences
iu_summaries_2002 = data.frame(IU_ID = haematobium_2002$IU_ID,
                               mean = apply(haematobium_2002[,2:ncol(haematobium_2002)], 1, mean, na.rm=T),
                               sd = apply(haematobium_2002[,2:ncol(haematobium_2002)], 1, sd, na.rm=T),
                               min = apply(haematobium_2002[,2:ncol(haematobium_2002)], 1, min, na.rm=T),
                               max = apply(haematobium_2002[,2:ncol(haematobium_2002)], 1, max, na.rm=T)) 

iu_summaries_2013 = data.frame(IU_ID = haematobium_2013$IU_ID,
                               mean = apply(haematobium_2013[,2:ncol(haematobium_2013)], 1, mean, na.rm=T),
                               sd = apply(haematobium_2013[,2:ncol(haematobium_2013)], 1, sd, na.rm=T),
                               min = apply(haematobium_2013[,2:ncol(haematobium_2013)], 1, min, na.rm=T),
                               max = apply(haematobium_2013[,2:ncol(haematobium_2013)], 1, max, na.rm=T)) 

iu_summaries_2022 = data.frame(IU_ID = haematobium_2022$IU_ID,
                               mean = apply(haematobium_2022[,2:ncol(haematobium_2022)], 1, mean, na.rm=T),
                               sd = apply(haematobium_2022[,2:ncol(haematobium_2022)], 1, sd, na.rm=T),
                               min = apply(haematobium_2022[,2:ncol(haematobium_2022)], 1, min, na.rm=T),
                               max = apply(haematobium_2022[,2:ncol(haematobium_2022)], 1, max, na.rm=T))

# high prev. IUs
high_prev_2002 = iu_summaries_2002 %>%
  filter((is.na(sd)))

high_prev_2013 = iu_summaries_2013 %>%
  filter((is.na(sd)))

high_prev_2022 = iu_summaries_2022 %>%
  filter((is.na(sd))) 

high_prev_ius = high_prev_2022 %>%
  filter(high_prev_2022$IU_ID %in% high_prev_2002$IU_ID & high_prev_2022$IU_ID %in% high_prev_2013$IU_ID)

# remove these IUs from maps
haematobium_2002 = haematobium_2002 %>%
  filter(!IU_ID %in% high_prev_ius$IU_ID)
haematobium_2013 = haematobium_2013 %>%
  filter(!IU_ID %in% high_prev_ius$IU_ID)
haematobium_2022 = haematobium_2022 %>%
  filter(!IU_ID %in% high_prev_ius$IU_ID)

############## MDA histories ############## 

# Read in MDA history
sch_histories_raw = read.csv(file.path(kPathToMapsSCH, "Schisto_IU_Cleaned_1.csv")) %>%
  mutate(EpiCov_binned = case_when(
    EpiCov <= 15 ~ 0,
    EpiCov > 15 & EpiCov <= 75 ~ 0.15,
    EpiCov > 75 ~ 0.75
  ))
sch_histories = sch_histories_raw %>%
  filter(IU_ID_MAPPING %in% haematobium_2002$IU_ID)  %>%
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
  filter(IU_ID_MAPPING %in% haematobium_2002$IU_ID)  %>%
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

add_burnin_year = data.frame(IU_ID_MAPPING = haematobium_2002$IU_ID,
                        Year = 1985,
                        TargetPop = "Adults",
                        Year_TargetPop = "1985_Adults",
                        EpiCov_binned = NA,
                        MDA_scheme = "Start burn in")
histories_joined = rbind(histories_joined,add_burnin_year) %>%
  arrange(IU_ID_MAPPING,Year)

histories_joined_full = data.frame(IU_ID_MAPPING = rep(haematobium_2002$IU_ID,each=41),
                              Year = rep(c(1985,rep(2003:2022,each=2)), length(haematobium_2002$IU_ID)),
                              TargetPop = rep(c("Adults",rep(c("SAC", "Adults"),20)),length(haematobium_2002$IU_ID))) %>%
  mutate(Year_TargetPop = paste0(Year,"_",TargetPop)) %>% 
  left_join(histories_joined) 
histories_joined_full$MDA_scheme[is.na(histories_joined_full$MDA_scheme)] = "Not delivered"


# Determine IUs with the same treatment history to assign batches

# find unique coverage histories
coverage_history = histories_joined_full %>%
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
all_treatments = data.frame(IU=haematobium_2002$IU_ID) %>%
  left_join(coverage_indexes) %>%
  left_join(MDA_scheme_indexes)

unique_all = all_treatments %>%
  group_by(across(c(-IU))) %>% 
  summarise(ius= list(IU)) %>%
  ungroup() %>%
  mutate(index = row_number())

# lookup table for IUs to batch IDs (used in original fitting)
iu_task_lookup = lapply(1:nrow(unique_all), function(i) data.frame(unique_all[["ius"]][[i]],unique_all[["index"]][i]))
iu_task_lookup = do.call(rbind,iu_task_lookup) 
colnames(iu_task_lookup) = c("IU_ID","TaskID")
iu_task_lookup = iu_task_lookup %>%
  left_join(unique(sch_histories %>% select(IU_ID_MAPPING)), by=c("IU_ID"="IU_ID_MAPPING"))

original_max_id = max(iu_task_lookup$TaskID)

# Overwrite MDA files for IUs with no treatments to add negligible coverage (otherwise python code will break)
# important to do this before reassigning batch numbers!
id_no_mda = original_max_id
iu_no_mda = unique((histories_joined_full %>% 
                      filter(IU_ID_MAPPING %in% (iu_task_lookup %>% 
                                                   filter(TaskID %in% id_no_mda))$IU_ID))$IU_ID_MAPPING)
inds = which(histories_joined_full$IU_ID_MAPPING %in% iu_no_mda & histories_joined_full$Year %in% c(2003,2004))
histories_joined_full[inds, colnames(histories_joined_full) %in% c("EpiCov_binned")] = 1e-11
histories_joined_full[inds, colnames(histories_joined_full) %in% c("MDA_scheme")] = notPZQ_name

# reassign problem IUs to their own batches for new fitting
# calculate group mean of tasks
batch_means_2002 = iu_task_lookup %>% 
  left_join(iu_summaries_2002, by=c("IU_ID")) %>%
  group_by(TaskID) %>% 
  summarise(mean_of_mean = mean(mean, na.rm=T))

batch_means_2013 = iu_task_lookup %>% 
  left_join(iu_summaries_2013, by=c("IU_ID")) %>%
  group_by(TaskID) %>% 
  summarise(mean_of_mean = mean(mean, na.rm=T))

batch_means_2022 = iu_task_lookup %>% 
  left_join(iu_summaries_2022, by=c("IU_ID")) %>%
  group_by(TaskID) %>% 
  summarise(mean_of_mean = mean(mean, na.rm=T))

# find IUs with low prevalences
# only keep low_prev_ius if max for that IU is less than mean of batch 
#threshold = 0.20
low_prev_2002 = iu_summaries_2002 %>%
  #filter(max < threshold) %>%
  left_join(iu_task_lookup, by="IU_ID") %>%
  left_join(batch_means_2002, by="TaskID") %>%
  filter(max < mean_of_mean)

low_prev_2013 = iu_summaries_2013 %>%
  #filter(max < threshold)  %>%
  left_join(iu_task_lookup, by="IU_ID") %>%
  left_join(batch_means_2013, by="TaskID") %>%
  filter(max < mean_of_mean)

low_prev_2022 = iu_summaries_2022 %>%
  #filter(max < threshold) %>%
  left_join(iu_task_lookup, by="IU_ID") %>%
  left_join(batch_means_2022, by="TaskID") %>%
  filter(max < mean_of_mean)

# find IUs and corresponding batches with low prev in any year
low_prev_ius = data.frame(IU_ID=sort(unique(c(low_prev_2022$IU_ID,low_prev_2013$IU_ID,low_prev_2002$IU_ID)))) %>%
  filter(!IU_ID %in% high_prev_ius$IU_ID) %>% 
  left_join(iu_task_lookup, by="IU_ID")
unique_low_prev_batches = data.frame(old_TaskID = sort(unique(low_prev_ius$TaskID))) %>%
  mutate(new_TaskID = max(iu_task_lookup$TaskID) + row_number())
# join new task ID to low_prev_ius
low_prev_ius = left_join(low_prev_ius, unique_low_prev_batches, by=c("TaskID" = "old_TaskID")) %>%
  arrange(TaskID)

# update batch IDs
iu_task_lookup$TaskID[iu_task_lookup$IU_ID %in% low_prev_ius$IU_ID] = low_prev_ius$new_TaskID
#View(iu_task_lookup %>% left_join(low_prev_ius,by="IU_ID") %>% filter(!TaskID.x == new_TaskID )) #sanity check: should have no rows

# save 2nd iteration of batches! 
save(iu_task_lookup, file=file.path(kPathToMapsArtefacts, "iu_task_lookup_haema.rds"))

# save batch IDs to maps file
haematobium_maps = list(list(data = iu_task_lookup %>% left_join(haematobium_2002, by=c("IU_ID"))),
                        list(data = iu_task_lookup %>% left_join(haematobium_2013, by=c("IU_ID"))),
                        list(data = iu_task_lookup %>% left_join(haematobium_2022, by=c("IU_ID"))))
save(haematobium_maps, file=file.path(kPathToMapsArtefacts, "haematobium_maps.rds"))

# Get into format for python model
coverage_wide = histories_joined_full %>%
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
drug_wide = unique(histories_joined_full %>%
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
kPathToEndgameInputs <- file.path(kPathToArtefacts, "endgame_inputs", "sch-haematobium")
if (!dir.exists(kPathToEndgameInputs)) {dir.create(kPathToEndgameInputs, recursive = TRUE)}

# Determine which batch IDs to process
batch_ids <- if (!is.null(opts$id)) {
  if (opts$id > max(iu_task_lookup$TaskID)) {
    stop(paste("Specified batch ID", opts$id, "exceeds maximum available batch ID", max(iu_task_lookup$TaskID)))
  }
  opts$id
} else {
  1:max(iu_task_lookup$TaskID)
}

print(paste("Processing batch IDs:", paste(batch_ids, collapse=", ")))

for (id in batch_ids){

  iu_file<-file.path(kPathToEndgameInputs,paste0("IUs_MTP_",id,".csv"))
  ius_per_batch = iu_task_lookup %>%
    filter(TaskID == id)
  ius = matrix(ius_per_batch$IU_ID,ncol=1)
  write.table(ius,file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
  
  for (iu in ius[1]) {
    
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
    mda_path<-file.path(kPathToEndgameInputs,paste0("InputMDA_MTP_",id,".xlsx"))
    write_xlsx(list("Platform Coverage" = mda_cov_iu_xl, "MarketShare" = mda_drug_iu_xl),path=mda_path, col_names = F) # write input MDA file
  }
}



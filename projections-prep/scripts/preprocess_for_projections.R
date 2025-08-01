# assuming working directory is post_AMIS_analysis

library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")

set.seed(100)

# uncomment for relevant species
#species <- "ascaris"
#failed_ids = c(1007,1127,1128,1333,252,269)
#failed_ids_sigma0.025 = c(1007,1127,1128,1333,252,269)

#species <- "hookworm"
#failed_ids = c(269,512)
#failed_ids_sigma0.025 = c(269,512)

species <- "trichuris"
failed_ids = c(1106,1168,1306,326,47,83)
failed_ids_sigma0.025 = c(1106,1168,1306,326,47,83)

#species <- "haematobium"
#failed_ids = c()
#failed_ids_sigma0.025 = c()

#species <- "mansoni_low_burden"
#failed_ids = c()
#failed_ids_sigma0.025 = c()

#species <- "mansoni_high_burden"
#failed_ids = c()
#failed_ids_sigma0.025 = c()

# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
if(species == "trichuris"){
  load("../Maps-STH/iu_task_lookup_trichuris.rds")
}else if (species %in% c("ascaris","hookworm")){
  load("../Maps-STH/iu_task_lookup_sth.rds")  
}else if (species == "haematobium"){
  load("../Maps-SCH/iu_task_lookup_haema.rds")
} else {
  load("../Maps-SCH/iu_task_lookup_mansoni.rds")
}
 
num_batches <- max(iu_task_lookup$TaskID)

cat(paste0("species: ",species, " \n"))
cat(paste0("num_batches: ",num_batches, " \n"))
  
ids_sample_pars = setdiff(1:num_batches,failed_ids_sigma0.025)


# this directory required for projections
if (!dir.exists("../ntd-model-sch/Man_MDA_vacc/")) {dir.create("../ntd-model-sch/Man_MDA_vacc/")}

#directory to save sampled parameters
InputPars_MTP_path_species <- paste0("InputPars_MTP_",species,"/")
if (!dir.exists(InputPars_MTP_path_species)) {dir.create(InputPars_MTP_path_species)}

# sample parameters and save draws
sampled_params_all = c()
for(id in ids_sample_pars){

  #### Load AMIS output (this is just to get IU names and initial ESS)
  if(!id %in% failed_ids){
    load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
  } else {
    load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
  }

  ess = amis_output$ess
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)

  if(!id %in% failed_ids){
    iu_names_lt200 = iu_names[ess<200]
    iu_names_ge200 = iu_names[!ess<200]
  } else {
    iu_names_lt200 = iu_names # if failed when sigma=0.0025 then use sigma=0.025 for all IUs
    iu_names_ge200 = NULL
  }

  #### Sample draws from the posterior
  num_sub_samples_posterior <- 200
	  
  for (iu in iu_names) {
    # use sigma=0.025 results if ESS < 200 when sigma=0.0025
    if(iu %in% iu_names_ge200 & (!id %in% failed_ids)){
      load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
      
    } else {
      load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
      
    }
    
    sampled_params <- sample_parameters(x = amis_output, n_samples = num_sub_samples_posterior, locations = which(iu_names==iu))

    # this puts all ius in the same folder
    file_name <- paste0(InputPars_MTP_path_species, paste0("InputPars_MTP_",iu,".csv"))
    write.csv(sampled_params, file=file_name, row.names = F)

    sampled_params_iu = cbind(IU_ID=iu,sampled_params)
    sampled_params_all = rbind(sampled_params_all,sampled_params_iu)

  }

  if(id%%100==0){cat(paste0("id=",id, "; "))}

}

save(sampled_params_all,file= paste0(InputPars_MTP_path_species,"InputPars_MTP_allIUs_",species,".rds"))
cat(paste0("Produced samples for all IUs in InputPars_MTP_",species, "/ \n"))


# Realocate files in correct format for projections

if(species %in% c("ascaris","hookworm","trichuris")){
  df_IU_country <- read.csv("../Maps-STH/table_iu_idx_STH.csv") # same for all species as jsut want country codes
} else if (species == "haematobium"){ 
  df_IU_country <- read.csv("../Maps-SCH/table_iu_idx_haematobium.csv")
} else {
  df_IU_country <- read.csv("../Maps-SCH/table_iu_idx_mansoni.csv")
}

countries <- sort(unique(df_IU_country$country))

#create directory
proj <- paste0("../ntd-model-sch/projections/")
if (!dir.exists(proj)) {dir.create(proj)}

path_species <- paste0("../ntd-model-sch/projections/",species,"/")
if (!dir.exists(path_species)) {dir.create(path_species)}

for(country in countries){
  path_country <- paste0("../ntd-model-sch/projections/",species,"/",country,"/")
  if (!dir.exists(path_country)) {dir.create(path_country)}
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
    path_iu <- paste0("../ntd-model-sch/projections/",species,"/",country,"/",country, iu0,"/")
    if (!dir.exists(path_iu)) {dir.create(path_iu)}
    
    file_name_old <- paste0(InputPars_MTP_path_species, paste0("InputPars_MTP_",iu,".csv"))
    sampled_params <- read.csv(file_name_old)
    
    sampled_params <- sampled_params[, c("seed", "R0", "k")]
    colnames(sampled_params) <- c("seed", "r0", "k")
    file_name_new <- paste0(path_iu, paste0("Input_Rk_",species_prefix, country, iu0,".csv"))
    write.csv(sampled_params, file=file_name_new, row.names = F)
    
  }
  
  if(id%%100==0){cat(paste0("id=",id, "; "))}

}

cat(paste0("Samples realocated for all IUs in InputPars_MTP_",species, "/ \n"))


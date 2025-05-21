# assuming working directory is post_AMIS_analysis

library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")


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

if(species %in% c("ascaris","hookworm","trichuris")){
  df_IU_country <- read.csv("../Maps-STH/table_iu_idx_STH.csv") # same for all species as jsut want country codes
} else if (species == "haematobium"){
  df_IU_country <- read.csv("../Maps-SCH/table_iu_idx_haematobium.csv")
} else {
  df_IU_country <- read.csv("../Maps-SCH/table_iu_idx_mansoni.csv")
}


# sample parameters and save draws
insufficient_ess = c()

for(id in ids_sample_pars){

  #### Load AMIS output
  if(file.exists(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata"))){
	  
  load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output

  ess = amis_output$ess
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  
  if(!id %in% failed_ids){
    iu_names_lt200 = iu_names[ess<200]
    iu_names_ge200 = iu_names[!ess<200]
  } else {
    iu_names_lt200 = iu_names # if failed when sigma=0.0025 then use sigma=0.025 for all IUs
    iu_names_ge200 = NULL
  }
  for(iu in iu_names_lt200){
    wh <- which(df_IU_country$IU_CODE==iu) 
    if(length(wh)!=1){stop("iu must be found exactly once in df_IU_country")} 
    country <- df_IU_country[wh, "country"] 
 
    iu0 <- paste0(country,sprintf("%05d", as.integer(iu)))
    insufficient_ess = c(insufficient_ess,iu0)
  }
 
  if(id%%100==0){cat(paste0("id=",id, "; "))}
  }
}

insufficient_ess_mat = data.frame(IU=insufficient_ess)
write.csv(insufficient_ess_mat, file=paste0("../ntd-model-sch/projections/",species,"/IUsWithInsufficientESS_",species,".csv"),row.names=F)

print("Produced IUsWithInsufficientESS.csv")



# assuming working directory is run-sch

library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")

## uncomment for relevant species
#species <- "ascaris"
#failed_ids = c(1007,1127,1128,1333,252,269)

#species <- "hookworm"
#failed_ids = c(269,512)

species <- "trichuris"
failed_ids = c(1106,1168,1306,326,47,83)

#species <- "haematobium"
#failed_ids = c()

#species <- "mansoni_low_burden"
#failed_ids = c()

#species <- "mansoni_high_burden"
#failed_ids = c()

# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
if(species == "trichuris"){
  load("../Maps-STH/iu_task_lookup_trichuris.rds")
}else if (species %in% c("ascaris","hookworm")){
  load("../Maps-STH/iu_task_lookup_sth.rds")  
}else if (species == "haematobium"){
  load("../Maps-SCH/iu_task_lookup_sch.rds")
  iu_task_lookup = iu_task_lookup_updated
  rm(iu_task_lookup_updated)
} else {
  load("../Maps-SCH/iu_task_lookup_mansoni.rds")
  iu_task_lookup = iu_task_lookup_mansoni
  rm(iu_task_lookup_mansoni)
}
 
num_batches <- max(iu_task_lookup$TaskID)
#num_batches <- max(iu_task_lookup$TaskID[which(iu_task_lookup$TaskID<10000)])

cat(paste0("species: ",species, " \n"))
cat(paste0("num_batches: ",num_batches, " \n"))
  
ids_sample_pars = setdiff(1:num_batches,failed_ids)

			  # save ESS 
ess_all_ius = c()
for(id in ids_sample_pars){
  
  load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  ess = amis_output$ess
  
  ess_all_ius = rbind(ess_all_ius, cbind(IU_ID = iu_names, TaskID=id, ess=ess))
}

ess_all_ius = as.data.frame(ess_all_ius) %>%
  mutate_if(is.character, as.numeric)

# find completed IUs that should use sigma=0.025 output
ess_ius_lt200 = ess_all_ius %>%
  filter(ess<200)

print("Batches with insufficient ESS:")
print(paste0(unique(ess_ius_lt200$TaskID),collapse=","))

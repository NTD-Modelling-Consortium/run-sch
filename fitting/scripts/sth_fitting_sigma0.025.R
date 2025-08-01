
# on cluster:
id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# for testing:
#id = 199
#setwd("~/Documents/STH-endgame-apr25/ntd-model-sch/")
#num_cores_to_use <- 8

# should be "ascaris", "hookworm" or "trichuris"
species = "trichuris"
task = "fitting"

library(dplyr)
library(AMISforInfectiousDiseases)

# assuming working directory is ntd-model-sch
source("../run-sch/R/amis_integration_sigma0.025.R")

args <- commandArgs(trailingOnly=TRUE)
num_cores_to_use <- parallel::detectCores()
if(length(args) == 1) {
    num_cores_to_use <- as.integer(args)
}

print(paste0('Using ', num_cores_to_use, ' cores'))

reticulate::use_virtualenv("../.venv/", required=TRUE)
sch_simulation <- get_amis_integration_package()
reticulate::py_config()

fixed_parameters <- sch_simulation$FixedParameters(
    # the higher the value of N, the more consistent the results will be
    # though the longer the simulation will take
    number_hosts = 500L,
    # no intervention
    coverage_file_name = ifelse(species=="trichuris",
                                paste0("endgame_inputs_STH/InputMDA_MTP_trichuris_",id,".xlsx"),
                                paste0("endgame_inputs_STH/InputMDA_MTP_",id,".xlsx")),
    demography_name = "UgandaRural",
    # cset the survey type to Kato Katz with duplicate slide
    survey_type = "KK2",
    parameter_file_name = paste0("STH_params/",species,"_params.txt"),
    coverage_text_file_storage_name = paste0("Man_MDA_vacc_",species,"_",id,"_sigma0.025.txt"),
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
# - Map samples are in discrete time, and we assume that samples for 2000 refer to the end of the year 2000 (approximately, 2000.99999 ~ 2001.0).
# - Therefore, if we have map samples for the (end of) years 2000, 2013, 2018, these should be compared to the simulations at times 2001.0, 2014.0, 2019.0 if we use the discrete-valued vector "year_indices".
# - Thus, we need to pass year_indices <- c(16L,29L,34L) to the Python code
# - The Python model will return a matrix with only 3 columns called 16L,29L,34L.
year_indices <- c(16L,29L,34L)

# Load prevalence map and filter rows for TaskID == id
load(paste0("../Maps-STH/",species,"_maps.rds"))
prevalence_map = get(paste0(species,"_map_allyears"))
prevalence_map = lapply(1:length(prevalence_map), function(t){
  output=list(data = as.matrix(prevalence_map[[t]]$data %>% 
                                 filter(TaskID==id) %>% 
                                 select(mu, sigma)), 
              likelihood = prevalence_map[[t]]$likelihood)
  rownames(output$data) = rownames(prevalence_map[[t]]$data)[prevalence_map[[t]]$data$TaskID==id]
  return(output)
  
})

# load prior
source(paste0("../run-sch/",species,"_prior.R"))
prior = Prior

# Algorithm parameters
amis_params<-default_amis_params()
amis_params$max_iters=50
amis_params$n_samples=1000
amis_params$target_ess=500
amis_params$sigma=0.025
amis_params$boundaries=c(-Inf,Inf)
amis_params$boundaries_param = matrix(c(R_lb,k_lb,R_ub,k_ub),ncol=2)

# shell to save trajectories
# commented out to pass Github tests (maybe this is bad practice...)
trajectories <- c() # save simulated trajectories as code is running
if (!dir.exists("../trajectories")) {dir.create("../trajectories")}
save(trajectories,file=paste0("../trajectories/trajectories_",id,"_",species,"_sigma0.025.Rdata"))

# Run AMIS
st<-Sys.time()
amis_output <- AMISforInfectiousDiseases::amis(
    prevalence_map,
    build_transmission_model(prevalence_map, fixed_parameters, year_indices, num_cores_to_use),
    prior,
    amis_params,
    seed = id
)
en<-Sys.time()
dur_amis<-as.numeric(difftime(en,st,units="mins"))
if (!dir.exists("../AMIS_output")) {dir.create("../AMIS_output")}
save(amis_output,file=paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata"))

print(amis_output)
summary(amis_output)
cat("--------------------- \n")
print(paste0("AMIS run time: ", round(dur_amis, digits=2), " minutes"))

# save summary
ess<-amis_output$ess
n_success<-length(which(ess>=amis_params[["target_ess"]]))
failures<-which(ess<amis_params[["target_ess"]])
n_failure<-length(failures)
if (n_failure>0) {cat(paste(failures,id,ess[failures]),file = paste0("../ESS_NOT_REACHED_",species,"_sigma0.025.txt"),sep = "\n", append = TRUE)}
if (!file.exists(paste0("../summary_",species,"_sigma0.025.csv"))) {cat("ID,n_failure,n_success,n_sim,min_ess,duration_amis,durarion_subsampling\n",file=paste0("../summary_",species,"_sigma0.025.csv"))}
cat(id,n_failure,n_success,length(amis_output$seeds),min(ess),dur_amis,NA,"\n",sep=",",file=paste0("../summary_",species,"_sigma0.025.csv"),append=TRUE)


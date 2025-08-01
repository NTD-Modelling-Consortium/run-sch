#!/bin/bash

#SBATCH --output running_hist_maps.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

#####################
# Your task is here #
#####################
# Add extra commands here to load a recent version of R
module purge
module load R/4.3.2-gfbf-2023a

stdbuf -i0 -o0 -e0 command
####################
# End of your task #
####################

# Now that you have loaded R above, we can run our R script
Rscript prepare_histories.R
Rscript prepare_histories_trichuris.R
Rscript prepare_maps_allspecies.R


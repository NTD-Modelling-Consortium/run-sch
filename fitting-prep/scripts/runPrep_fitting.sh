#!/bin/bash

#SBATCH --output log/outputs-fitting-prep.out-%A_%a
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

cd ./Maps-SCH

# Load modules
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript prepare_histories_and_maps_haematobium.R 
Rscript prepare_histories_and_maps_mansoni.R 



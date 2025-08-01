#!/bin/bash

#SBATCH --output log/find-low-ess.out-%A
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

# Load modules
module purge
module load R/4.3.2-gfbf-2023a

stdbuf -i0 -o0 -e0 command

# Change directory 
cd run-sch 

# Run R script
Rscript find_lowESS_ids.R


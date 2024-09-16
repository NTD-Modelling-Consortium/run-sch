#!/bin/bash

#SBATCH --output log/mtp-sch-haematobium.out-%A_%a
#SBATCH --array=1-1126
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00

# Change directory
cd run-sch

# Load modules
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1
source .venv/bin/activate
unset RETICULATE_PYTHON

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript sch_fitting.R 12


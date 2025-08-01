#!/bin/bash

#SBATCH --output log/mtp-sth-trichuris-sigma0.025.out-%A_%a
#SBATCH --array=1106,1168,1306,326,47,83,1365,1155
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00

# Change directory
cd ntd-model-sch

# Load modules
module purge
module load R/4.3.2-gfbf-2023a

source ../.venv/bin/activate
unset RETICULATE_PYTHON

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript ../run-sch/sth_fitting_sigma0.025.R 12


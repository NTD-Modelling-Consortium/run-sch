#!/bin/bash

#SBATCH --output log/mtp-sth-hookworm.out-%A_%a
#SBATCH --array=1-1387
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00

#1-1387

# Change directory
cd ntd-model-sch

# Load modules
module purge
module load R/4.3.2-gfbf-2023a

source ../.venv/bin/activate
unset RETICULATE_PYTHON

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript ../run-sch/sth_fitting.R 12


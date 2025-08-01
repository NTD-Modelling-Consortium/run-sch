#!/bin/bash

#SBATCH --output log/mtp-sth-ascaris-sigma0.025.out-%A_%a
#SBATCH --array=1333,1248,1317,1366,1387,749,1056,1106,1139,1220
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00

#1007,1127,1128,252,269,5,40,55,85,90,94,106,113,117,118,141,160,171,190,206,209,225,267,299,335,346,387,397,409,449,451,454,484,485,509,522,523,771
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


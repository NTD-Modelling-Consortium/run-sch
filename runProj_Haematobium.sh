#!/bin/bash

#SBATCH --output log/proj_haematobium.out-%A_%a
#SBATCH --array=
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=6:00:00

#5248 total

module purge
module load Python/3.11.3-GCCcore-12.3.0

stdbuf -i0 -o0 -e0 command

cd ntd-model-sch/
source ../.venv/bin/activate

python ../run-sch/sth_amis/sch_projections_per_IU.py ${SLURM_ARRAY_TASK_ID}


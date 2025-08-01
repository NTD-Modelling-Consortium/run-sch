#!/bin/bash

#SBATCH --output log/proj_hookworm.out-%A_%a
#SBATCH --array=0-5359
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=6:00:00

#index starts from 0-5359

module purge
module load Python/3.11.3-GCCcore-12.3.0

stdbuf -i0 -o0 -e0 command

cd ntd-model-sch/
source ../.venv/bin/activate

python ../run-sch/sth_amis/sth_projections_per_IU.py ${SLURM_ARRAY_TASK_ID}


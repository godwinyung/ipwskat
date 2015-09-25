#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH -p serial_requeue
#SBATCH --mem=4000
#SBATCH -o output/SLURMout/sim%a.out
#SBATCH -e output/SLURMerror/sim%a.err

R CMD BATCH --quiet --no-restore --no-save "--args taskID=$SLURM_ARRAY_TASK_ID interval=10" sim-setting.r output/Rout/sim${SLURM_ARRAY_TASK_ID}.Rout
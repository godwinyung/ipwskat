#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-00:10
#SBATCH -p serial_requeue
#SBATCH --mem=200

python ../code/genhap.py out.hap-1 regions3e3.txt haplotypes3e3 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID
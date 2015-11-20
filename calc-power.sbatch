#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH -p serial_requeue
#SBATCH --mem=400
#SBATCH -o output/SLURMout/sim%a.out
#SBATCH -e output/SLURMerror/sim%a.err

python ../../code/calcpower.py results/summary/t1e*.pvalues
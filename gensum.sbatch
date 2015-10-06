#!/bin/bash
#
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-00:10
#SBATCH -p serial_requeue
#SBATCH --mem=200
#
for i in $(seq 1 2); do
  for j in $(seq 1 2); do
    python ../../code/gensum.py t1e${i}${j} results/summary/t1e${i}${j}
  done
done
#!/bin/sh

#SBATCH -J bin
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH -p general
#SBATCH -o genhap.out
#SBATCH -e genhap.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=godwin.y.yung@gmail.com
#SBATCH --mem=600

python genhap.py

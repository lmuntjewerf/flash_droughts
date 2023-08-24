#!/bin/bash
# 
# SLURM batch script 
#

#SBATCH --job-name=calc_FDs
#SBATCH --qos=nf
#SBATCH --time=24:00:00
#SBATCH --mem=16000            # memory per node in MB 

python calc_FDs.py
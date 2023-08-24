#!/bin/bash
#
# SLURM batch script 
#

#SBATCH --job-name=Rnet_PET
#SBATCH --qos=nf
#SBATCH --time=24:00:00
#SBATCH --mem=16000            # memory per node in MB 

#-------------------------------
# commands to be executed
#-------------------------------

scenario=$1
member=$2
frequency=$3

cd /home/nklm/Px_drought/PET

./run_Rnet_PET_calc_args.sh $scenario $member $frequency




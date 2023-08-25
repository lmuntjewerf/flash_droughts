#!/bin/bash
#
# SLURM batch script 
#

#SBATCH --job-name=PET
#SBATCH --qos=nf
#SBATCH --time=24:00:00
#SBATCH --mem=16000            # memory per node in MB 

#-------------------------------
# commands to be executed
#-------------------------------

# scenario=$1
# member=$2
# frequency=$3

cd /home/nklm/Px_drought/flashdroughts/PET_PM/

# ./run_PET_calc_args.sh $scenario $member $frequency

./run_PET_calc_args.sh 


#!/bin/bash

#set -ex

#-------------------
# calculate PET for daily ERA5 output
# this script takes arguments:
#
#
#-------------------

diri='/scratch/nklm/Px_flashdroughts/PET_data/'
startyear=2020
endyear=2022

for year in $(seq $startyear $endyear); do 
  python3 calc_PET.py --diri $diri -y $year
done
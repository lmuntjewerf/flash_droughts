#!/bin/bash

#set -ex

#-------------------
# calculate PET for daily ERA5 output
# this script takes arguments:
# diri
# diro 
# year
#-------------------

diri='/scratch/nkkw/Karin/P2_flashdroughts/meteodata_ERA5/gridded/'
diro='/scratch/nklm/Px_flashdroughts/PET_ERA5/'
startyear=1940
endyear=2022

for year in $(seq $startyear $endyear); do 
  python3 calc_PET.py --diri $diri --diro $diro -y $year
done
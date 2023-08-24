#!/bin/bash

#set -ex

#-------------------
# calculate PET for daily ERA5 output
# this script takes arguments:
#
#
#-------------------

diri=
freq='day'
startyear=2020
endyear=2022

for year in $(seq $startyear $endyear); do 
  python3 calc_Rnet_PET.py --diri $diri --freq $freq -y $year
done
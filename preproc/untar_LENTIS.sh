#!/bin/bash

set -ex

#-------------------
# LENTIS data that is downloaded from ECFS. 
# this script first move the tarzip files into the right subfolders, before untarring them. 
#
# untar LENTIS data in folder
# delete tar.gz files 
#-------------------


# freq='Amon'
# declare -a vars
# vars=(pet pr)

freq='Lmon'
declare -a vars
vars=(mrsos mrso)

# freq='day'
# declare -a vars
# vars=(mrsos mrso)
# # vars=(pet pr)

# switches for scenario to untar [do (1) or don't (0)]
scenario_h=0
scenario_s=1

# -------- no user edits below this license

# this is the folder where the tar files are: 
LENTISSCRATCH=/perm/nklm/Px_drought/testing/
# this is the folder where the untarred files go: 
FILOTOGO=/perm/nklm/Px_drought/testing/LENTIS/ec-earth/cmorised_by_var/


# the depth of the tar files is variable, due to how I saved it. 
# to get the depth, we need to list the contents of a tar or tar.gz file without extracting it and grep a netcdf file so we can count the depth of the folder structure: 
# tar -tvf my-data.tar.gz 'search-pattern'
# file_to_extract=ec/res4/scratch/nklm/cmorisation/cmorised-results/EC-EARTH-AOGCM/h010/CMIP6/CMIP/KNMI/EC-Earth3/historical/r1i1p5f1//Amon/pet/gr/v20220601/pet_mon_EC-Earth3_historical_r1i1p5f1_gr_200101-200112.nc
# depth=$(awk -F/ '{print NF-1}' <<< "$file_to_extract")


# depth=18

# -------- historical
for var in ${vars[@]}; do
 if [ $var == 'pet' ]
 then
   depth=18
 else
   depth=3
 fi
 echo $depth

 if (( $scenario_h ))
 then
    echo; echo " *** UNTARRING hxxx ***"; echo
    scenario=h
    directory=${FILOTOGO}${scenario}xxx/${freq}/${var}/
    mkdir -p ${directory}
    cd ${directory}
    for i in {1..1}; do
      for j in {0..0}; do
         mkdir -p ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
         tar -xvzf ${LENTISSCRATCH}/${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var}.tar.gz --strip-components=$depth -C  ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
         chmod -R u+w ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
        #  rm ${LENTISSCRATCH}/${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var}.tar.gz
      done ;
    done
 fi


 # -------- scenario_future
 if (( $scenario_s ))
 then
  echo; echo " *** UNTARRING sxxx ***"; echo
  scenario=s
  directory=${FILOTOGO}${scenario}xxx/${freq}/${var}/
  mkdir -p ${directory}
  cd ${directory}
  for i in {1..1}; do
    for j in {0..0}; do
       mkdir -p ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
       tar -xvzf ${LENTISSCRATCH}/${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var}.tar.gz --strip-components=$depth -C  ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
       chmod -R u+w ${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var} ;
    #    rm ${LENTISSCRATCH}/${scenario}"$(printf "%02d" $i)"${j}_${freq}_${var}.tar.gz
    done ;
  done
 fi


done
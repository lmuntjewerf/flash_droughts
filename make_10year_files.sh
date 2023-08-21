#!/bin/bash

set -ex

#-------------------
# postprocess LENTIS data into files of 10 years 
# adjusts lon axis 
#
# puts in KNMI-MSO-varEX structure in the workstation
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



# switches for scenario to postprocess [do (1) or don't (0)]
scenario_h=1
scenario_s=0


# -------- no user edits below this line

# this is the folder where the tar files are: 
LENTISSCRATCH=/perm/nklm/Px_drought/testing/LENTIS/
# this is the folder where the untarred files go: 


TEMPDIR=/scratch/nklm/temp/
mkdir -p ${TEMPDIR}



if (( $scenario_h ))
 then
   scenario=h
   CMIP_scen=historical
   twoletter=PD
#    CMOR_NO=1
   startyear=2000
   endyear=2009
fi

if (( $scenario_s ))
 then
   scenario=s
   CMIP_scen=ssp245
   twoletter=2K
#    CMOR_NO=5
   startyear=2075
   endyear=2084
fi

for var in ${vars[@]}; do
  ORIG_DIR=${LENTISSCRATCH}/ec-earth/cmorised_by_var/${scenario}xxx/${freq}/${var}/
  OPRUIM_DIR=${LENTISSCRATCH}/${twoletter}/${freq}/${var}/
  mkdir -p ${OPRUIM_DIR}

  if [ $var == 'pet' ] && [ $freq == 'Amon' ]
  then
   filefreq='mon'
  else
   filefreq=$freq
  fi
  echo $filefreq

  for P in {1..1}; do
    for E in {0..0}; do
      CODE=${scenario}"$(printf "%02d" $P)"${E}

      member_dir=${ORIG_DIR}/${CODE}_${freq}_${var}/

    #   # SET CMOR CONVENTION
    #   CMOR=r${CMOR_NO}"$(printf "%02d" $P)"${E}i${E}p5f1
    #   if (( $scenario_h )) && (( ${P} == "1" )) && (( ${E} == "0" ));
    #   then
    #           CMOR=r1i1p5f1   # note first member different convention
    #   fi

    #   echo ${CMOR}


      if [ ${freq} = "day" ] || [ ${freq} = "3hr" ] # If you want to say OR use double pipe (||).
      then
              years=($(seq ${startyear} ${endyear}))
              for year in ${years[@]}; do
                      echo ${year}
                      cdo sellonlatbox,-180,180,-90,90 ${member_dir}/${var}_${filefreq}_EC-Earth3_${CMIP_scen}_*_gr_${year}0101-${year}1231.nc ${TEMPDIR}/${CODE}_${var}_${freq}_${year}.nc
              done
              cdo mergetime ${TEMPDIR}/${CODE}_${var}_${freq}_*.nc ${OPRUIM_DIR}/${var}_${CODE}.nc
              rm -rf ${TEMPDIR}/${CODE}_${var}_${freq}_*.nc
      elif [ ${freq} = "fx" ] || [ ${freq} = "Ofx" ] # If you want to say OR use double pipe (||).
      then
              cdo sellonlatbox,-180,180,-90,90 ${member_dir}/${var}_${filefreq}_EC-Earth3_${CMIP_scen}_*_*.nc ${OPRUIM_DIR}/${var}.nc
      else
              echo 'nothing'
              cdo mergetime ${member_dir}/${var}_${filefreq}_EC-Earth3_${CMIP_scen}_*_??_*.nc ${TEMPDIR}/${var}_${CODE}.nc
              cdo sellonlatbox,-180,180,-90,90 ${TEMPDIR}/${var}_${CODE}.nc ${OPRUIM_DIR}/${var}_${CODE}.nc
              rm -rf ${TEMPDIR}/${var}_${CODE}.nc
      fi


    done ;
  done ;
done
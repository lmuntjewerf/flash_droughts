#!/bin/bash


# ik heb nodig: 
# - tmean
# - tmin
# - tmax
# - wind
# - rh --> maken in python
# - rs 

# units --> maken in python
# - temperatuur in degC
# - rsds in MJ/day
# - rh in %
# - wind in m/s

# 2 dit wil ik voor de jaren 2020-2022


datadir_in='/scratch/nkkw/ERA5/'
datadir_out='/scratch/nklm/Px_flashdroughts/PET_data/test_2D/'

# for wind
# cdo -chname,ws,sfcWind -mergetime ${datadir_in}/raw3_202?.nc ${datadir_out}/sfcWind.nc

# for tmin tmax dewpointtemp
# cdo mergetime ${datadir_in}/raw2_202?.nc ${datadir_out}/raw2.nc
# cdo -chname,mx2t,tmax -selvar,mx2t ${datadir_out}/raw2.nc ${datadir_out}/tmax.nc
# cdo -chname,d2m,tdew -selvar,d2m ${datadir_out}/raw2.nc ${datadir_out}/tdew.nc
# cdo -chname,mn2t,tmin -selvar,mn2t ${datadir_out}/raw2.nc ${datadir_out}/tmin.nc
# rm ${datadir_out}/raw2.nc
# for i in  0 1 2; do 
#   cdo -chname,mx2t,tasmax -selvar,mx2t ${datadir_in}/raw2_202${i}.nc ${datadir_out}/tasmax_202${i}.nc
#   cdo -chname,d2m,tdew -selvar,d2m ${datadir_in}/raw2_202${i}.nc ${datadir_out}/tdew_202${i}.nc
#   cdo -chname,mn2t,tasmin -selvar,mn2t ${datadir_in}/raw2_202${i}.nc ${datadir_out}/tasmin_202${i}.nc
# done 
# cdo  mergetime ${datadir_out}/tasmin_202${i}.nc ${datadir_out}/tasmin.nc
# cdo  mergetime ${datadir_out}/tasmax_202${i}.nc ${datadir_out}/tasmax.nc
# cdo  mergetime ${datadir_out}/tdew_202${i}.nc ${datadir_out}/tdew.nc

for var in  tas tasmin tasmax  tdew sfcWind rsds ; do
  cdo ydaymean ${datadir_out}/${var}.nc ${datadir_out}/${var}_ydaymean.nc
done

# for i in  0 1 2; do 
#   cdo -chname,ssrd,rsds -selvar,ssrd ${datadir_in}/raw_202${i}.nc ${datadir_out}/rsds_202${i}.nc
#   cdo -chname,t2m,tas -selvar,t2m ${datadir_in}/raw_202${i}.nc ${datadir_out}/tas_202${i}.nc
# done 
# cdo  mergetime ${datadir_out}/rsds_202${i}.nc ${datadir_out}/rsds.nc
# cdo  mergetime ${datadir_out}/tas_202${i}.nc ${datadir_out}/tas.nc

# for var in tas tmin tmax sfcWind rsds tdew; do
#   cdo daymean ${datadir_out}/${var}.nc ${datadir_out}/${var}_daymean.nc
# done


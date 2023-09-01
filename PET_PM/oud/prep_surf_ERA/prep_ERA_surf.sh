#!/bin/bash

# bash script that prepares 2D gridded lat and elevation from ERA5 for PET PM claculations
# input: geopotential on 0.1 res grid. 
# steps: 
# 1. devide geopotential with gravitational constant to get surface elevation
# 2. write out lat as a 2D field


PETprep_dir='/perm/nklm/Px_flashdroughts/PET/prep/'

# make ERA5 griddes and write to txt file
infile='/scratch/nkkw/ERA5/raw_1940.nc'
cdo -griddes ${infile} >> ${PETprep_dir}ERAgridfile.txt

# regrid the surf_para.nc file 
surf_para_file='surf_para.nc'
surf_para_file_regridded='surf_para_regridded.nc'
cdo -remapbil,${PETprep_dir}ERAgridfile.txt ${PETprep_dir}${surf_para_file} ${PETprep_dir}${surf_para_file_regridded}
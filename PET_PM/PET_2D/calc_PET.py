#!/usr/bin/env python
# coding: utf-8
import xarray as xr
import pyet
import glob
import os
import argparse
import sys

import numpy as np

"""
#              CALC PET from ERA5 reanalysis data
# --------------------------------------------------------------------------
# Script to calculate potential evapotranspiration following the ASCE Penman-Monteith  method
# as offered by the pyet package https://github.com/pyet-org/pyet/
# This script writes the output to netcdf while maintaining the file attributes. Output in diri. 
#
# Requires the following cmorised data:
# rsds,tas,tasmin,tasmax,sfcWind,hurs,orog
# 
#
# Author: Laura Muntjewerf (laura.muntjewerf@knmi.nl)
#
# Call as follow to calculate PET:
#
#      ./calc_PET.py --diri 'input dir' --year 'year'
#
"""




def calc_pet_tonetcdf(diri,diro,year):
    """
    Calculate Potential evapotranspiration with ASCE Penman-Monteith from cmorised data. 
    With the package pyet: Estimation of Potential Evapotranspiration [https://github.com/pyet-org/pyet]
    This function assumes cmorised file structure
    
    In: path up before the frequency to tas, tasmax, tasmin, sfcWind, hurs, elevation. Assumes cmorised file structure (diri example: /scratch/nklm/cmorisation/cmorised-results/EC-EARTH-AOGCM/ft30/CMIP6/CMIP/KNMI/EC-Earth3/historical/r1i1p5f1/day/ )
    In: year to calculate pet for
    
    Out: netcdf file saved in f'{diri}/{freq}/pet (daily, monthly, gridded or timeseries) 
    """
    import xarray as xr
    import pyet

    import numpy as np
    
    import os
    import glob
    
    # ==============================
    # Do calculations


    static_file = '/perm/nklm/Px_flashdroughts/PET/prep/surf_para_regridded.nc'
    gravitational_constant = 9.81


    elevation      = xr.open_dataset(static_file)['z'][0,:,:]/gravitational_constant
    lat1d          = xr.open_dataset(static_file)['latitude']* np.pi / 180
    lat2d          = lat1d.expand_dims(dim={"longitude":elevation.longitude}, axis=1) 
    del lat1d

    tmean       = xr.open_dataset(f'{diri}/tas_{year}.nc/',engine='netcdf4')["tas"]
    tasmin      = xr.open_dataset(f'{diri}/tasmin_{year}.nc',engine='netcdf4')["tasmin"]
    tasmax      = xr.open_dataset(f'{diri}/tasmax_{year}.nc',engine='netcdf4')["tasmax"]
    wind        = xr.open_dataset(f'{diri}/sfcWind_{year}.nc',engine='netcdf4')["sfcWind"]
    rs          = xr.open_dataset(f'{diri}/rsds_{year}.nc',engine='netcdf4')["rsds"] * (86400 / 1000000) #/(60*60)
    hurs        =  xr.open_dataset(f'{diri}/hurs_{year}.nc',engine='netcdf4')["hurs"]


    # tmean =xr.open_mfdataset(f'{diri}/{freq}/tas/*/*/tas*{year}*.nc')["tas"]-273.15
    # tmin  =xr.open_mfdataset(f'{diri}/{freq}/tasmin/*/*/tasmin*{year}*.nc')["tasmin"]-273.15
    # tmax  =xr.open_mfdataset(f'{diri}/{freq}/tasmax/*/*/tasmax*{year}*.nc')["tasmax"]-273.15
    # wind  =xr.open_mfdataset(f'{diri}/{freq}/sfcWind/*/*/sfcWind*{year}*.nc')["sfcWind"]
    # rh    =xr.open_mfdataset(f'{diri}/{freq}/hurs/*/*/hurs*{year}*.nc')["hurs"]
    # #rn    =xr.open_mfdataset(f'{diri}/{freq}/rnet/*/*/rnet*{year}*.nc')["rnet"] * 86400 / 1000000  # concert to [MJ/m2day]
    # rs    =xr.open_mfdataset(f'{diri}/{freq}/rsds/*/*/rsds*{year}*.nc')["rsds"] * 86400 / 1000000  # concert to [MJ/m2day] 

    # tmean =xr.open_dataset(f'{diri}/tas_Rhine.nc/',engine='netcdf4')["tas"]
    # tmin  =xr.open_dataset(f'{diri}/tasmin_Rhine.nc',engine='netcdf4')["tasmin"]
    # tmax  =xr.open_dataset(f'{diri}/tasmax_Rhine.nc',engine='netcdf4')["tasmax"]
    # wind  =xr.open_dataset(f'{diri}/sfcWind_Rhine.nc',engine='netcdf4')["sfcWind"]
    # rh    =xr.open_dataset(f'{diri}/hurs_Rhine.nc',engine='netcdf4')["hurs"]
    # rs    =xr.open_dataset(f'{diri}/rsds_Rhine.nc',engine='netcdf4')["rsds"] * (86400 / 1000000)  # concert to [MJ/m2 day] 

    # elevation    = 477 #xr.open_mfdataset(f'{diri}/fx/orog/gr/v20220601/orog_fx_EC-Earth3_*_gr.nc')["orog"]
    # lat          = 50 * np.pi / 180 # xr.open_mfdataset(f'{diri}/fx/orog/gr/v20220601/orog_fx_EC-Earth3_*_gr.nc')["lat"]* np.pi / 180
    # # lat2         = lat.expand_dims(dim={"lon":elevation.lon}, axis=1)
    
    # tmin['lat']=tmean['lat']
    # tmax['lat']=tmean['lat']
    # wind['lat']=tmean['lat']
    # rh['lat']=tmean['lat']
    # rs['lat']=tmean['lat']
    # elevation['lat']=tmean['lat']

    # Do calculations
    pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs, elevation=elevation, lat=lat2d, tmax=tasmax, tmin=tasmin, rh=hurs.data)
   
    del tmean
    del tasmin
    del tasmax
    del wind
    del hurs
    del rs
    del elevation
    del lat2d
    
    # ==============================
    # Get/make the right folder and file name to save the PET data in
    filename=f'{diri}/tas_{year}.nc/'
    print(filename)
    ds=xr.open_dataset(filename,engine='netcdf4')

    time_attrs=ds.time.attrs
    # coords={'time': (['time'], ds.time.data,time_attrs)}
    lat_attrs=ds.latitude.attrs
    lon_attrs=ds.longitude.attrs
    coords={'time': (['time'], ds.time.data,time_attrs),
            'latitude': (['latitude'], ds.latitude.data,lat_attrs),
            'longitude': (['longitude'], ds.longitude.data,lon_attrs),
            }


        
    pet_filename=f'pet_{year}.nc'
        
    
    # ==============================
    # Prepare xarray.Dataset to save as Netcdf file

    # define global attributes
    attrs=ds.attrs
    
   
    # define data with variable attributes
    data_vars = {'pet_pm_asce':(['time', 'latitude', 'longitude'],  pet_pm_asce.data, 
                         {'units': 'mm day-1', 
                          'standard_name':'pet_ASCE_PM',
                          'long_name':'Potential evapotranspiration with ASCE Penman-Monteith',
                          'comment':'The surface potential evapotranspiration with the ASCE Penman-Montheith method is calculated with the package pyet: Estimation of Potential Evapotranspiration [https://github.com/pyet-org/pyet] using the following ERA5 reanalysis variables: daily tmean, daily wind speed, incoming solar radiation in [MJ/m2day], daily tmax, daily tmin, daily rh, fixed elevation and fixed latitude in the following equation: pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)'})}

    # data_vars = {'pet_pm_asce':(['time'],  pet_pm_asce.data, 
    #                      {'units': 'mm day-1', 
    #                       'standard_name':'pet_ASCE_PM',
    #                       'long_name':'Potential evapotranspiration with ASCE Penman-Monteith',
    #                       'comment':'The surface potential evapotranspiration with the ASCE Penman-Montheith method is calculated with the package pyet: Estimation of Potential Evapotranspiration [https://github.com/pyet-org/pyet] using the following ERA5 reanalysis variables: daily tmean, daily wind speed, incoming solar radiation in [MJ/m2day], daily tmax, daily tmin, daily rh, fixed elevation and fixed latitude in the following equation: pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)'})}



    # create dataset
    ds_new = xr.Dataset(data_vars=data_vars, 
                     coords=coords, 
                     attrs=attrs)
    
    # ==============================
    # Save rnet in netcdf in right location
    # mind time encoding, otherwise it will throw error when opening with nvciew --> (ncview: netcdf_dim_value: unknown data type (10) for dimension time)
    ds_new.to_netcdf(diro+pet_filename, unlimited_dims='time', encoding={'time': {'dtype': 'i4'}})
    


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate PET for postprocessed ERA5 files for specified year.')
    parser.add_argument('-y', '--year', help='Year to do the calculations for.')
    parser.add_argument('-i', '--diri', help='Input directory')
    parser.add_argument('-o', '--diro', help='Output directory')
    
    args = parser.parse_args()

    if args.year != None and args.diri != None:
        diri=args.diri
        diro=args.diro
        year=args.year
        print("Start PET calc for year "+str(args.year)+" from "+str(args.diri)+" & saving file in"+str(args.diro))
        calc_pet_tonetcdf(diri,diro,year)
        

if __name__ == '__main__':
    sys.exit(main())


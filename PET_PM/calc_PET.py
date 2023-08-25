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
#              CALC Rnet and PET for cmorised model data
# --------------------------------------------------------------------------
# Script to calculate net surface radiation and potential evapotranspiration 
# writing these to netcdf while maintaining the file attributes. Output in diri. 
#
# Requires the following cmorised data:
# rsds,rsus,rlds,rlus,tas,tasmin,tasmax,sfcWind,hurs,orog
# 
#
# Author: Laura Muntjewerf (laura.muntjewerf@knmi.nl)
#
# Call as follow to calculate Rnet and PET:
#
#      ./perturb_var.py --diri 'input dir' --freq 'cmor freq; day or mon' --year 'year'
#
"""




def calc_pet_tonetcdf(diri, year):
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

    # tmean =xr.open_mfdataset(f'{diri}/{freq}/tas/*/*/tas*{year}*.nc')["tas"]-273.15
    # tmin  =xr.open_mfdataset(f'{diri}/{freq}/tasmin/*/*/tasmin*{year}*.nc')["tasmin"]-273.15
    # tmax  =xr.open_mfdataset(f'{diri}/{freq}/tasmax/*/*/tasmax*{year}*.nc')["tasmax"]-273.15
    # wind  =xr.open_mfdataset(f'{diri}/{freq}/sfcWind/*/*/sfcWind*{year}*.nc')["sfcWind"]
    # rh    =xr.open_mfdataset(f'{diri}/{freq}/hurs/*/*/hurs*{year}*.nc')["hurs"]
    # #rn    =xr.open_mfdataset(f'{diri}/{freq}/rnet/*/*/rnet*{year}*.nc')["rnet"] * 86400 / 1000000  # concert to [MJ/m2day]
    # rs    =xr.open_mfdataset(f'{diri}/{freq}/rsds/*/*/rsds*{year}*.nc')["rsds"] * 86400 / 1000000  # concert to [MJ/m2day] 

    tmean =xr.open_dataset(f'{diri}/tas_Rhine.nc/',engine='netcdf4')["tas"]
    tmin  =xr.open_dataset(f'{diri}/tasmin_Rhine.nc',engine='netcdf4')["tasmin"]
    tmax  =xr.open_dataset(f'{diri}/tasmax_Rhine.nc',engine='netcdf4')["tasmax"]
    wind  =xr.open_dataset(f'{diri}/sfcWind_Rhine.nc',engine='netcdf4')["sfcWind"]
    rh    =xr.open_dataset(f'{diri}/hurs_Rhine.nc',engine='netcdf4')["hurs"]
    rs    =xr.open_dataset(f'{diri}/rsds_Rhine.nc',engine='netcdf4')["rsds"] * (86400 / 1000000)  # concert to [MJ/m2 day] 

    elevation    = 477 #xr.open_mfdataset(f'{diri}/fx/orog/gr/v20220601/orog_fx_EC-Earth3_*_gr.nc')["orog"]
    lat          = 50 * np.pi / 180 # xr.open_mfdataset(f'{diri}/fx/orog/gr/v20220601/orog_fx_EC-Earth3_*_gr.nc')["lat"]* np.pi / 180
    # lat2         = lat.expand_dims(dim={"lon":elevation.lon}, axis=1)
    
    # tmin['lat']=tmean['lat']
    # tmax['lat']=tmean['lat']
    # wind['lat']=tmean['lat']
    # rh['lat']=tmean['lat']
    # rs['lat']=tmean['lat']
    # elevation['lat']=tmean['lat']

    

    pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs.values, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh.data)
   
    del tmean
    del tmin
    del tmax
    del wind
    del rh
    del rs
    del elevation
    del lat
    
    # ==============================
    # Get/make the right folder and file name to save the rnet data in
    filename=glob.glob(f'{diri}/tas_Rhine.nc/')
    ds=xr.open_dataset(filename)

    time_attrs=ds.time.attrs
    coords={'time': (['time'], ds.time.data,time_attrs)}
    # lat_attrs=ds.lat.attrs
    # lon_attrs=ds.lon.attrs
    # coords={'time': (['time'], ds.time.data,time_attrs),
    #         'lat': (['lat'], ds.lat.data,lat_attrs),
    #         'lon': (['lon'], ds.lon.data,lon_attrs)
    #         }

    
    # filename=glob.glob(f'{diri}/{freq}/tas/*/*/tas*{year}*.nc')
    path=filename[0].split('tas')
    
    pet_path=path[0]+'pet'+path[1]
    isExist = os.path.exists(pet_path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(pet_path)
        
    pet_filename='pet'+path[-1]
        
    
    # ==============================
    # Prepare xarray.Dataset to save as Netcdf file

    # define global attributes
    attrs=ds.attrs
    
   
    # define data with variable attributes
    data_vars = {'pet_pm_asce':(['time', 'lat', 'lon'],  pet_pm_asce.data, 
                         {'units': 'mm day-1', 
                          'standard_name':'pet_ASCE_PM',
                          'long_name':'Potential evapotranspiration with ASCE Penman-Monteith',
                          'comment':'The surface potential evapotranspiration with the ASCE Penman-Montheith method is calculated with the package pyet: Estimation of Potential Evapotranspiration [https://github.com/pyet-org/pyet] using the following LENTIS variables: daily tmean (tas_day), daily wind speed (sfcWind_day), incoming solar radiation in [MJ/m2day] (rs = rsds_day * 86400 / 1000000), daily tmax (tasmax_day), daily tmin (tasmin_day), daily rh (hurs_day) fixed elevation (orog_fx) and fixed latitude (ds.lat) in the following equation: pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)'})}

    data_vars = {'pet_pm_asce':(['time'],  pet_pm_asce.data, 
                         {'units': 'mm day-1', 
                          'standard_name':'pet_ASCE_PM',
                          'long_name':'Potential evapotranspiration with ASCE Penman-Monteith',
                          'comment':'The surface potential evapotranspiration with the ASCE Penman-Montheith method is calculated with the package pyet: Estimation of Potential Evapotranspiration [https://github.com/pyet-org/pyet] using the following LENTIS variables: daily tmean (tas_day), daily wind speed (sfcWind_day), incoming solar radiation in [MJ/m2day] (rs = rsds_day * 86400 / 1000000), daily tmax (tasmax_day), daily tmin (tasmin_day), daily rh (hurs_day) fixed elevation (orog_fx) and fixed latitude (ds.lat) in the following equation: pet_pm_asce = pyet.pm_asce(tmean, wind, rs=rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)'})}



    # create dataset
    ds_new = xr.Dataset(data_vars=data_vars, 
                     coords=coords, 
                     attrs=attrs)
    
    # ==============================
    # Save rnet in netcdf in right location
    # mind time encoding, otherwise it will throw error when opening with nvciew --> (ncview: netcdf_dim_value: unknown data type (10) for dimension time)
    ds_new.to_netcdf(pet_path+pet_filename, unlimited_dims='time', encoding={'time': {'dtype': 'i4'}})
    


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate PET and Rnet for files in cmorized diri for specified year.')
    parser.add_argument('-y', '--year', help='Year to do the calculations for.')
    parser.add_argument('-d', '--diri', help='Input directory cmor structure up')
    
    args = parser.parse_args()

    if args.year != None and args.diri != None:
        diri=args.diri
        year=args.year
        #print("Start Rnet calc for year "+str(args.year)+" in "+str(args.diri)+" ...")
        print("Start PET calc for year "+str(args.year)+" in "+str(args.diri)+" ...")
        calc_pet_tonetcdf(diri,year)
        

if __name__ == '__main__':
    sys.exit(main())


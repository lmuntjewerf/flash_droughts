#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import xarray as xr
import sys

from datetime import datetime, timedelta

# def get_info_from_filename(fili):
#     '''
#     Purpose: string splitting to extract some information from the input file name.

#     Input: File name of CSV file with all variations in it (string)
#         example: 'Rhine_ESI14_dwindow21_jump-2.0.csv'
#     Output: (list of strings, len=4)
#         basin: River basin (example: Rhine)
#         indexscale: index + scale (example: ESI14)
#         window: amount of days over which we take the jump (example: 21)
#         jump: jump in index (example: 2.0)
#     '''
#     basin, indexscale = fili.split('_')[0:-2]
#     window = int(fili.split('_')[-2][-2::])
#     jump = int(fili.split('_')[-1].split('.')[0][-2::])/10
#     return basin, indexscale, window, jump

def get_info_from_filename(fili):
    '''
    Purpose: string splitting to extract some information from the input file name.

    Input: File name of CSV file with all variations in it (string)
        example: 'Rhine_ESI14_dwindow21_jump-2.0.csv'
    Output: (list of strings, len=4)
        basin: River basin (example: Rhine)
        indexscale: index + scale (example: ESI14)
        window: amount of days over which we take the jump (example: 21)
        jump: jump in index (example: 2.0)
        end_threshold: end threshold of index (example: -1.5)
    '''
    basin, indexscale = fili.split('_')[0:-3]
    window = int(fili.split('_')[-3][1::])
    jump = int(fili.split('_')[-2].split('.')[0][-2::])/10
    end_threshold = int(fili.split('_')[-1].split('.')[0][-2::])/10
    return basin, indexscale, window, jump, end_threshold

def check_leap_year(year):
    leap_year = False
    if((year % 400 == 0) or  (year % 100 != 0) and  (year % 4 == 0)):  
        leap_year = True
    return leap_year


def define_ts_length(year,month,day,window,extra_len):

    # make the right type
    window = int(window)
    extra_len = int(extra_len)
    
    start_FD_time = datetime(year,month,day)

    # when do we want the composite timeseries to start and end:
    start_ts_time = datetime(year,month,day) + timedelta(days=-(extra_len))
    end_ts_time = datetime(year,month,day) + timedelta(days=(window+extra_len))

    if check_leap_year(year+1): 
        leap_day = datetime(year+1,2,29)
        # is next years leap day between the start of the timeseries and the start of the event: add a day BEFORE
        print(start_ts_time)
        print(leap_day)
        print(start_FD_time)
        print(end_ts_time)
        if start_ts_time <= leap_day <= start_FD_time:
            # this is always false, in principe, because we are in the year BEFORE the leapday. 
            start_ts_time = datetime(year,month,day) + timedelta(days=-(extra_len+1))
        elif start_FD_time <= leap_day <= end_ts_time:
            # this can be true for any start_FD_time day between 01-12 and 31-12 in the year BEFORE the leapday.  
            end_ts_time = datetime(year,month,day) + timedelta(days=(window+extra_len+1))

    elif check_leap_year(year):
        leap_day = datetime(year,2,29)

        # is the leap day between the start of the timeseries and the start of the event: add a day BEFORE
        if start_ts_time <= leap_day <= start_FD_time:
            start_ts_time = datetime(year,month,day) + timedelta(days=-(extra_len+1))
        # is the leap day between the the start of the event and the end of the timeseries: add a day AFTER
        elif start_FD_time <= leap_day <= end_ts_time:
            end_ts_time = datetime(year,month,day) + timedelta(days=(window+extra_len+1))



    # # however, if there is a leap day in this period, we need to add an extra day to the time series before or after  
    # if check_leap_year(year):
    #     leap_day = datetime(year,2,29)
        
    #     # is the leap day between the start of the timeseries and the start of the event: add a day BEFORE
    #     if start_ts_time <= leap_day <= start_FD_time:
    #         start_ts_time = datetime(year,month,day) + timedelta(days=-(extra_len+1))
    #     # is the leap day between the the start of the event and the end of the timeseries: add a day AFTER
    #     elif start_FD_time <= leap_day <= end_ts_time:
    #         end_ts_time = datetime(year,month,day) + timedelta(days=(window+extra_len+1))

    
    # # TODO verbeter leap year voor kerst case 
    # # diri_csv = '/perm/nklm/Px_flashdroughts/ERA5_FD_events/'
    # # fili_csv = 'Rhine_SMI-7_w28_j20_e15.csv'
    # # event 33 



    return start_ts_time, end_ts_time


# def make_composites(df, window, extra_len, indexscale, basin):
#     '''
    
#     '''
#     # make the right type
#     window = int(window)
#     extra_len = int(extra_len)

#     # initialise timeseries
#     timeseries = []

#     for index, row in df.iterrows():

#         # step 1: gather what info we have about the start date of the FD events:
#         year=row["year"]
#         month = row["month"]
#         day = row["day"]

#         # step 2: compute the dates we want the composite timeseries to start and end
#         start_ts_time, end_ts_time = define_ts_length(year, month, day, window, extra_len)

#         # step 3: open the data file and extract the timeseries we want for the event we want
#         diri = '/scratch/nklm/Px_flashdroughts/indices_ERA5/'
#         fili = f'{indexscale}_{basin}.nc'
#         ds=xr.open_dataset(diri+fili)

#         # da = ds['ESI14'].sel(time=slice(start_ts_time,end_ts_time))
#         da = ds[indexscale].sel(time=slice(start_ts_time,end_ts_time))
#         dt = pd.to_timedelta(range(-(extra_len),int(window)+extra_len), unit='D')

#         try: 
#             da = da.assign_coords(time = dt.days)
#             da = da.assign_coords(event = index)
#             timeseries.append(da)
#         except ValueError:
#             pass


#     ds_combined = xr.concat(timeseries, dim='event')
#     return ds_combined

def make_composites(df, var, window, extra_len, indexscale, basin, anom=False):
    '''
    
    '''
    # make the right type
    window = int(window)
    extra_len = int(extra_len)

    # initialise timeseries
    timeseries = []

    for index, row in df.iterrows():

        # step 1: gather what info we have about the start date of the FD events:
        year=row["year"]
        month = row["month"]
        day = row["day"]

        # step 2: compute the dates we want the composite timeseries to start and end
        start_ts_time, end_ts_time = define_ts_length(year, month, day, window, extra_len)

        # step 3: open the data file and extract the timeseries we want for the event we want
        indexscales = []
        # even op een rijtje zetten elke we hebben
        for scale in (7,14,21,28):
            for FDindex in ['SPI','SPEI','ESI','SMI']:
                indexscale = f'{FDindex}-{scale}'
                indexscales.append(indexscale)
        if var in indexscales:
            diri = '/scratch/nklm/Px_flashdroughts/indices_ERA5/'
            fili = f'{var}_{basin}.nc'
            ds=xr.open_dataset(diri+fili)
        elif var in ['pr','et','pet','tas','rsds','mrsos']:
            # diri = '/scratch/nkkw/Karin/P2_flashdroughts/meteodata_ERA5/'
            # if anom: 
            #     filis = f'{var}_{basin}_????_anom.nc'
            #     ds=xr.open_mfdataset(diri+filis)
            # else: 
            #     filis = f'{var}_{basin}_????.nc'
            #     ds=xr.open_mfdataset(diri+filis)
            diri = '/scratch/nklm/Px_flashdroughts/meteodata_ERA5/'
            if anom: 
                fili = f'{var}_{basin}_anom.nc'
                ds=xr.open_dataset(diri+fili)
            else: 
                fili = f'{var}_{basin}.nc'
                ds=xr.open_dataset(diri+fili)

        # da = ds['ESI14'].sel(time=slice(start_ts_time,end_ts_time))
        da = ds[var].sel(time=slice(start_ts_time,end_ts_time))
        dt = pd.to_timedelta(range(-(extra_len),int(window)+extra_len), unit='D')

        try: 
            da = da.assign_coords(time = dt.days)
            da = da.assign_coords(event = index)
            timeseries.append(da)
        except ValueError:
            pass


    ds_combined = xr.concat(timeseries, dim='event')
    return ds_combined

def save_as_netcdf(ds, var, basin, index, scale, window, jump, diro,end_threshold, anom=False):
    '''
    Purpose: Save ds of composite data to netcdf in right location for var == varname
    mind encoding, otherwise it will throw error when opening with nvciew 
         --> (ncview: netcdf_dim_value: unknown data type (10) for dimension time) - this is an integer (days)
         --> (ncview: ncview: netcdf_dim_value: unknown data type (12) for dimension ensemble_member_year) - this is a string of 9 characters (h010_2000 etc.)
    
    Input: ds and strings to compose the output filename and location

    Out: netcdf file saved in f'{diro}/{filo}' 
    '''
    jumpout = int(abs(jump)*10)
    end_thresholdout = int(abs(end_threshold)*10)

    if anom: 
        filo = f'{var}_anom_{basin}_{index}-{scale}_w{window}_j{jumpout}_e{end_thresholdout}.nc'
    else: 
        filo = f'{var}_{basin}_{index}-{scale}_w{window}_j{jumpout}_e{end_thresholdout}.nc'

    ds.to_netcdf(diro+filo, unlimited_dims='event', encoding={'time': {'dtype': 'i4'}, 'event': {'dtype': 'i4'}})


def main():
    
    diri = '/perm/nklm/Px_flashdroughts/ERA5_FD_events/'
    diro = '/scratch/nklm/Px_flashdroughts/composites/'
    basin = 'Rhine'
    extra_len = int(65)
    end_threshold = -1.5

    # vars = ['SPI-14', 'tas']
    # vars = ['ESI-14',]
    # vars = ['pr','et','pet','tas','rsds','mrsos']
    # vars = ['SMI-7', 'SPEI-28', 'SPI-28']
    vars = ['ESI-14',]
    anom = False


    #for index in ['SPI','SPEI','ESI','SMI']:
    for index in ['SPEI','ESI',]:
        # for scale in [7,14,21,28]:
        for scale in [14,21]:
            # for window in [7,14,21,28]:
            for window in [21,28]:
                for jump in np.arange(-1.5,-2.5,-0.5): 
                    jumpout = int(abs(jump)*10)
                    end_thresholdout = int(abs(end_threshold)*10)
                    fili = f'{basin}_{index}-{scale}_w{window}_j{jumpout}_e{end_thresholdout}.csv'
                    df = pd.read_csv(diri+fili,index_col=0)
                    basin, indexscale, window, jump, end_threshold = get_info_from_filename(fili)
                    # basin, indexscale, window, jump = get_info_from_filename(fili)

                    for var in vars: 
                        ds = make_composites(df, var, window, extra_len, indexscale, basin, anom)
                        save_as_netcdf(ds, var, basin, index, scale, window, jump, diro, end_threshold, anom)
                    


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python
# coding: utf-8


import xarray as xr
import pandas as pd
import numpy as np
import sys


# ## Goal: define flash drougths
# 
# In this notebook, I am working to a script that defines the flash drought from an indexes and a number of criteria. 
# 
# 
# ### The criteria we need (input): 
# 
# - we choose an index **index** of a certain time scale **scale**
# - we compute the first order derivate of the index with a rolling-mean window **window** 
# - the starting index value has to be higher than **start_threshold**
# - the end index value has to be lower than **end_threshold**
# - the first order derivative of the index must experience a jump in value of minimal size **jump**
# 
# 
# ### Output of the script:
# 
# - CSV file with start dates of flash droughts
# - file name must contain the values of the above input criteria

# ### steps: 
#  1. read in Index netcdf
#  2. calculate first order derivative
#  3. find where data fits the following conditions:
#     1. index(t) > start_threshold
#     2. index(t+window) < end_threshold
#     3. d_index(t+window) < jump (because the timeseries start at the time when the running window ends)

def read_in_ERA5_index(index, diri, basin, scale):
    # ds = xr.open_dataset(f'{diri}/ESI7_test.nc{index}{scale}_{basin}.nc')
    ds = xr.open_dataset(f'{diri}/{index}{scale}_test.nc')
    da = ds[f'{index}{scale}']
    return da


def first_order_deriv_rolling_sum(da, window):
    dda = da[0:-1].copy(deep=True)
    dda.values = pd.Series(np.diff(da)).rolling(window=window).sum()
    
    return dda



def compute_FDs(da, dda, window, start_threshold, end_threshold, jump):
    '''
    da is the index timeseries
    dda is the time derivate of da, with rolling window
    window is the rolling window of the dda
    start_threshold is below with value da should start at the beginning of the flash drought
    end_threshold is below which value da should end after #window days after beginning of flash drought
    jump is the minimal value of dda at #window days after beginning of flash drought
    '''
    flash_drought_dates = []
    for i in da.time.data:
        i_end = i+np.timedelta64(window,'D')
        if da.sel(time=i,method='nearest') > start_threshold and da.sel(time=i_end,method='nearest') < end_threshold and dda.sel(time=i_end,method='nearest') < jump:
            flash_drought_dates.append(i.astype('datetime64[D]'))
    df = pd.DataFrame(data={'FD_startdate':flash_drought_dates})
    return df 



def save_df_tocsv(df, diro, basin, index, scale, window, jump):
    filo = f'{basin}_{index}{scale}_dwindow{window}_jump{jump}.csv'
    print(diro+filo)
    df.to_csv(diro+filo)  


def main():
    diri = '/scratch/nklm/Px_flashdroughts/indices_ERA5/'
    diro = '/perm/nklm/Px_flashdroughts/ERA5_FD_events/'
    basin = 'Rhine'
    start_threshold = 0
    end_threshold = -1
    for index in ['SPI','SPEI','ESI','SMI']:
        for scale in [7,14,21,28]:
            for window in [7,14,21,28]:
                for jump in np.arange(-1.5,-5,-0.5): 
                    da = read_in_ERA5_index(index, diri, basin, scale)
                    dda = first_order_deriv_rolling_sum(da, window)
                    df = compute_FDs(da, dda, window, start_threshold, end_threshold, jump)
                    save_df_tocsv(df, diro, basin, index, scale, window, jump)



if __name__ == '__main__':
    sys.exit(main())


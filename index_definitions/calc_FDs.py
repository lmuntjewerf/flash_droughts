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
    ds = xr.open_dataset(f'{diri}/{index}-{scale}_{basin}.nc')
    da = ds[f'{index}-{scale}']
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

    # cluster events that are consecutive days
    df['event'] = df['FD_startdate'].diff().dt.days.ne(1).cumsum()
    # group these events in a new df
    df_temp = df.groupby('event')['FD_startdate'].agg(['min', 'max'])

    # take individual events that are window+7 days apart from each other and cluster them
    df_temp['indiv_event'] = df_temp['min'].diff().dt.days.ge(window+7).cumsum()

    # group these events in a new df
    df_new = df_temp.groupby('indiv_event')['min'].agg(['min'])
    df_new['FD_startdate'] = df_new['min']
    df_new = df_new.drop(columns=['min'])
    df_new['year'] = df_new.FD_startdate.map(lambda x: x.year)
    df_new['month'] = df_new.FD_startdate.map(lambda x: x.month)
    df_new['day'] = df_new.FD_startdate.map(lambda x: x.day)
    
    return df_new 



def save_df_tocsv(df, diro, basin, index, scale, window, jump):
    jumpout = int(abs(jump)*10)
    filo = f'{basin}_{index}-{scale}_w{window}_j{jumpout}.csv'
    print(diro+filo)
    df.to_csv(diro+filo)  


def main():
    diri = '/scratch/nklm/Px_flashdroughts/indices_ERA5/'
    diro = '/perm/nklm/Px_flashdroughts/ERA5_FD_events/'
    basin = 'Rhine'
    start_threshold = 0
    end_threshold = -1
    # for index in ['SPI','SPEI','ESI','SMI']:
    for index in ['ESI',]:
        # for scale in [7,14,21,28]:
        for scale in [14,21]:
            # for window in [7,14,21,28]:
            for window in [21,28]:
                for jump in np.arange(-1.5,-2.5,-0.5): 
                    da = read_in_ERA5_index(index, diri, basin, scale)
                    dda = first_order_deriv_rolling_sum(da, window)
                    df = compute_FDs(da, dda, window, start_threshold, end_threshold, jump)
                    save_df_tocsv(df, diro, basin, index, scale, window, jump)



if __name__ == '__main__':
    sys.exit(main())


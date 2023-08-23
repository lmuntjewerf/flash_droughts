#!/usr/bin/env python
# coding: utf-8

# ## Goal: Indexes for the flash drougth project
# 
# In this notebook, I am working to a script the calculates indexes we need for the flash drougth project. These will be timeseries for the area-integrated Rhine basin. 
# 
# 
# ### The indexes we need are: 
# 
# - SPI standardized precipitation index. Variable: precipitation
# - SPEI standardized precipitation/evapotranspiration index. Variable: precipitation - potential evapotranspiration
# - ESI/SESR evaporative stress index/standardized evaporative stress ratio. Variable: evapotranspiration/potential evapotranspiration
# - SMI soil moisture index. Variable: soil moisture
# 
# ### Input variables for the script
# 
# - Variable to standardize 
# - Time scale (For now we test all variables for 7, 14, 21, 28 day time scales)

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
# import functions (and packages) from R
r_time_series = robjects.r('ts')

def read_in_ERA5(var, diri, basin):
    ds = xr.open_mfdataset(f'{diri}/{var}_{basin}_????.nc')
    da = ds[var]
    return da


def make_index_timeseries(var, diri,basin):
    if var == 'pr' or var == 'mrsos':
        da = read_in_ERA5(var, diri, basin)
    elif var == 'wb':
        da_pr = read_in_ERA5('pr', diri, basin)
        da_pet = read_in_ERA5('pet', diri, basin)
        da = da_pr - da_pet
    elif var == 'es':
        da_et = read_in_ERA5('et', diri, basin)
        da_pet = read_in_ERA5('pet', diri, basin)
        da_et[da_et < 0 ] = 0
        da_pet[da_pet < 0] = 0
        da = da_et / da_pet
        da[da > 5] = 5        
    return da


def calc_standardized_index_daily(da, scale):
    SPEI_package = importr('SPEI')
    r_spei_function = robjects.r['spei']

    # convert to R timeseries
    r_da = r_time_series(robjects.FloatVector(da.values), start = robjects.IntVector([da.time.dt.year[0].values, da.time.dt.month[0].values]), frequency = 365)
    r_standardized_index = r_spei_function(r_da, scale=scale, na_rm=True, distribution='log-Logistic',verbose=False)


    # put standardized index in DataArray
    da_standardized_index = da.copy(deep=True)
    da_standardized_index.values = xr.DataArray(pandas2ri.ri2py_vector(r_standardized_index.rx2('fitted')))
    
    return da_standardized_index


def calc_index_to_netcdf(var, diri, diro, basin, scale):

    da = make_index_timeseries(var, diri,basin)
    da_standardized_index = calc_standardized_index_daily(da, scale)

    # define data with variable attributes as Xarray dataset

    if var == 'pr':
        indexname = f'SPI-{scale}'
        explanation = 'precipitation'
    elif var == 'mrsos':
        indexname = f'SMI-{scale}'
        explanation = 'top layer soil moisture'
    elif var == 'wb':
        indexname = f'SPEI-{scale}'
        explanation = 'water balance (pr - PET)'
    elif var == 'es':
        indexname = f'ESI-{scale}'
        explanation = 'evaporative stress (ET - PET)'


    var_attr = {'units': '-', 
                        'standard_name': f'{indexname}_{scale}d',
                        'long_name': f'{indexname} standardized index of {var} ({explanation}) on scale: {scale}d. ',
                        }

    coords={'time': (['time'], da_standardized_index.time.data, da_standardized_index.time.attrs)}
              
    ds_new = xr.Dataset(
    data_vars=dict(
        index=(['time'],  
                      da_standardized_index.data, 
                       var_attr,
    )),
    coords=coords
    )

    # da_standardized_index.name = indexname
    # da_standardized_index.attrs = var_attr
    # ds_new = da_standardized_index.to_dataset

    ds_new.time.encoding = {'zlib': False,
                        'shuffle': False,
                        'complevel': 0,
                        'fletcher32': False,
                        'contiguous': False,
                        'dtype': np.dtype('float64'),
                        'units': 'days since 1850-01-01 00:00:00',
                        'calendar': 'proleptic_gregorian'
                        }

    ds_new.index.encoding = {'_FillValue': np.nan,
                'missing_value':np.nan,
                'dtype': np.dtype('float64'),
                'zlib': False
            }

    # ds_new['index'].rename(indexname)
    ds_new[indexname] = ds_new['index']
    ds_new = ds_new.drop(['index'])

    ds_new.to_netcdf(f'{diro}/{indexname}_{basin}.nc')



def main():
    basin = 'Rhine'
    diri='/scratch/nkkw/Karin/P2_flashdroughts/meteodata_ERA5/'
    diro = '/scratch/nklm/Px_flashdroughts/indices_ERA5'
    for var in ['pr','mrsos','wb']:
        # for var in ['es',]:
        for scale in [7,14,21,28]:
            #for scale in [7,28]:
            calc_index_to_netcdf(var, diri, diro, basin, scale)


if __name__ == '__main__':
    sys.exit(main())

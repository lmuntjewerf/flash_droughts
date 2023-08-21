#!/usr/bin/env python
# coding: utf-8


import xarray as xr
import numpy as np

import sys


def open_one_LENTIS(var, freq, timeslice, i, j, diri, plev=None):
    """
    Open one LENTIS data file from the ensemle
    """
    if timeslice == 'PD': 
        letter = 'h'
    elif timeslice == '2K': 
        letter = 's'
        
    ens_member=f'{letter}{str(i).zfill(2)}{str(j)}'
    
    file=f"{diri}/{timeslice}/{freq}/{var}/{var}_{ens_member}.nc"
    
    if plev != None:
        ds=xr.open_dataset(file).sel(plev=plev)
    else:
        ds=xr.open_dataset(file)
    
    return ds, ens_member


def calc_boxstat(ds, var, area='global'): 
    """
    Compute spatial weighted mean
    ds      :  xarray DataArray
    """ 
    box_seasons = ds.time.dt.season

    if hasattr(ds, 'lat'):
            weights = np.cos(ds.lat * np.pi / 180)
    elif hasattr(ds, 'latitude'):
            weights = np.cos(ds.latitude * np.pi / 180)

    if area=='global':
        if hasattr(ds, 'lat'):
            boxstat = ds[var].weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].weighted(weights).mean(dim=('latitude','longitude')) 
    elif area=='nh':
        if hasattr(ds, 'lat'): 
            boxstat = ds[var].sel(lat=slice(0,90)).weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].sel(latitude=slice(0,90)).weighted(weights).mean(dim=('latitude','longitude')) 
    elif area=='europe':
        if hasattr(ds, 'lat'):
            boxstat = ds[var].sel(lat=slice(30,70),lon=slice(-10,40)).weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].sel(latitude=slice(30,70),longitude=slice(-10,40)).weighted(weights).mean(dim=('latitude','longitude'))  

    # regions for Pieter Slomp
    elif area in ['rhine','donau','euti','ganges','yangtze','nile','chad','niger','orinoco','mississippi']:
        path_mask = '/nobackup/users/slomp/basin_masks/' + area + '.nc'
        basin_mask = xr.open_dataset(path_mask).catchmentID[0,:,:]

        if hasattr(ds, 'lat'):
            boxstat = ds[var].where(basin_mask).weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].where(basin_mask).weighted(weights).mean(dim=('latitude','longitude'))   
            
    # regions for Leonie Hemelrijk
    elif area=='LH_westEU':
        if hasattr(ds, 'lat'):
            boxstat = ds[var].sel(lat=slice(30,60),lon=slice(-30,30)).weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].sel(latitude=slice(30,60),longitude=slice(-30,30)).weighted(weights).mean(dim=('latitude','longitude'))    
    elif area=='LH_rhine':
        if hasattr(ds, 'lat'):
            boxstat = ds[var].sel(lat=slice(47,55),lon=slice(2,10)).weighted(weights).mean(dim=('lat','lon')) 
        elif hasattr(ds, 'latitude'):
            boxstat = ds[var].sel(latitude=slice(47,55),longitude=slice(2,10)).weighted(weights).mean(dim=('latitude','longitude'))

    else:
        raise ValueError('unknown area: '+area)
            
            
    return boxstat, box_seasons

def calc_seasonal_avg(boxstat, box_seasons, seas='ANN'):
    '''
    '''
    if seas == 'ANN':
        seasonal_avg = boxstat.groupby('time.year').mean('time')
    elif seas in ['DJF', 'MAM','JJA','SON']:
        seasonal_avg = boxstat.sel(time=box_seasons==seas).groupby('time.year').mean("time")
    else:
        raise ValueError('unknown seas: '+seas)
        
    return seasonal_avg


def get_attr_info(var, freq, timeslice, diri):
    """
    
    """
    ds, ens_member = open_one_LENTIS(var, freq, timeslice, 1, 1, diri)
    
    # get global attributes
    attrs=ds.attrs
    
    # get coordinates
    time_attrs=ds.time.attrs
    lat_attrs=ds.lat.attrs
    lon_attrs=ds.lon.attrs

    # get coordinate bounds 
    time_bnds=ds.time_bnds
    lat_bnds=ds.lat_bnds
    lon_bnds=ds.lon_bnds
    
    # get part of the var attributes
    var_units=ds[var].units
    var_standard_name=ds[var].standard_name
    var_long_name=ds[var].long_name
    
    
    return attrs, time_attrs, time_bnds, var_units, var_standard_name, var_long_name


def get_avg_timebnds(var, freq, timeslice, diri, seas='ANN'):
    '''
    '''
    ds, ens_member = open_one_LENTIS(var, freq, timeslice, 1, 1, diri)
    
    if seas == 'ANN':
        time_bnds_str = ds.time_bnds.sel(bnds=0).groupby('time.year').min(dim='time')
        time_bnds_end = ds.time_bnds.sel(bnds=1).groupby('time.year').max(dim='time')
        time_avg = ds.time.groupby('time.year').mean('time')
    elif seas == 'DJF':
        time_bnds_str = ds.time_bnds.sel(bnds=0,time=ds.time.dt.season=="DJF").groupby('time.year').max("time")
        time_bnds_end = ds.time_bnds.sel(bnds=1,time=ds.time.dt.season=="DJF").groupby('time.year').min("time")
        time_avg = ds.time.sel(time=ds.time.dt.season=="DJF").groupby('time.year').mean("time")
    elif seas in [ 'MAM','JJA','SON']:
        time_bnds_str = ds.time_bnds.sel(bnds=0,time=ds.time.dt.season==seas).groupby('time.year').min("time")
        time_bnds_end = ds.time_bnds.sel(bnds=1,time=ds.time.dt.season==seas).groupby('time.year').max("time")
        time_avg = ds.time.sel(time=ds.time.dt.season==seas).groupby('time.year').mean("time")
    else:
        raise ValueError('unknown seas: '+seas)
        
    time_bnds_combined = np.vstack((time_bnds_str, time_bnds_end)).T
    
    return time_bnds_combined, time_avg


def calc_ens_statistics_and_to_netcdf(var, freq, timeslice, plev, diri, diro, area='global', seas='ANN', statistic='mean'):
    '''
    """Compute the annual-averaged, global weighted mean 
    for a given variable, for all ensemble members of a time slice"""
    '''
    
    # define variables of the computed quantity and its ensemble member
    ens_stat_values=[]
    ens_member_list=[]
    
    # Do the calculation
    for i in np.arange(1,16+1):
        for j in np.arange(0,9+1):
            ds, ens_member = open_one_LENTIS(var, freq, timeslice, i, j, diri, plev)
            #ds, ens_member = open_one_LENTIS(var, freq, timeslice, i, j, diri)
            boxstat,box_seasons = calc_boxstat(ds, var, area)
            seasonal_avg = calc_seasonal_avg(boxstat,box_seasons, seas)

            ens_stat_values.append(seasonal_avg)
            ens_member_list.append(ens_member)
            del ds,boxstat,seasonal_avg,ens_member

    
    # Prepare xarray.dataset to save as Netcdf file
    attrs, time_attrs, time_bnds, var_units, var_standard_name, var_long_name = get_attr_info(var,freq,timeslice, diri)
    
    # time_attrs = {'standard_name': 'time',
    #                 'long_name': 'time',
    #                 'bounds': 'time_bnds',
    #                 'axis': 'T'}
    
    time_bnds_combined, time_avg = get_avg_timebnds(var, freq, timeslice, diri, seas)
    
    ens_attrs = {'standard_name': 'ens_mem',
                 'long_name': 'Ensemble member',
                 'comment': 'Postprocessed by Laura Muntjewerf (KNMI). A time slice of KNMI-LENTIS consists of 16 ensemble member of 10 years. All simulations have a unique ensemble member label that reflects the forcing, and how the initial conditions are generated. The initial conditions have two aspects: the parent simulation from which the run is branched (macro perturbation, there are 16), and the seed relating to a particular micro-perturbation in the initial three-dimensional atmosphere temperature field (there are 10). The ensemble member label thus is a combination of: forcing (h for present-dsy/historical and s for +2K/SSP2-4.5), parent ID (number between 1 and 16), micro perturbation ID (number between 0 and 9) '}
    coords={'time': (['time'], time_avg.data, time_attrs),
            'ens': (['ens'], ens_member_list, ens_attrs)
            }
    
    var_attr = {'units': var_units, 
                        'standard_name': f'{var_standard_name}_{area}_{seas}_{statistic}',
                        'long_name': f'{statistic} {var_long_name} of {area} ({seas})'}
                    
    
    # define data with variable attributes as Xarray dataset
    ds_new = xr.Dataset(
    data_vars=dict(
        time_bnds=(["time", "bnds"], 
                     time_bnds_combined,{}),
        ens_bnds=(["ens"], 
                     ens_member_list,{}),
        var=(['ens','time'],  
                      ens_stat_values, 
                       var_attr,
    )),
    coords=coords,
    attrs=attrs,
    )
    
    # assign encoding
    # otherwise it will throw error when opening with nvciew --> (ncview: netcdf_dim_value: unknown data type (10) for dimension time) en (unknown data type (12) for dimension ens)
    ds_new.ens.encoding = {'zlib': False,
                     'shuffle': False,
                     'complevel': 0,
                     'fletcher32': False,
                     'contiguous': False,
                     'chunksizes': None,
                     'dtype': '|S1'}

    ds_new.time.encoding = {'zlib': False,
                        'shuffle': False,
                        'complevel': 0,
                        'fletcher32': False,
                        'contiguous': False,
                        'dtype': np.dtype('float64'),
                        'units': 'days since 1850-01-01 00:00:00',
                        'calendar': 'proleptic_gregorian'
                        }

    

    # create dataset
    # Save in right location
    
    if plev != None:
        filo=f'{timeslice}_ensemble_{area}_{seas}_{freq}_{var}_{plev}_mean.nc'
    else:
        filo=f'{timeslice}_ensemble_{area}_{seas}_{freq}_{var}_mean.nc'
    

    ds_new.to_netcdf(f'{diro}/{filo}', unlimited_dims='ens')
    

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate PET and Rnet for files in cmorized diri for specified year and freq.')
    parser.add_argument('-v', '--variable', help='Variable to do the calculations for.')
    parser.add_argument('-f', '--frequency', help='Freq to do the calculations for (Amon, Lmon).')
    parser.add_argument('-t', '--timeslice', help='Timeslice to do the calculations for (PD, 2K)')
    parser.add_argument('-i', '--diri', help='Input directory where the data is.')
    parser.add_argument('-o', '--diro', help='Output directory where to save the data')
    parser.add_argument('-a', '--area', help='Area to aggregate over (default = global)(options: global, nh, europe)', default='global', nargs='?')
    parser.add_argument('-s', '--seas', help='Season to calculte (default = ANN)(options: ANN, DJF, MAM, JJA, SON)', default='ANN', nargs='?')
    parser.add_argument('-c', '--stat', help='Statistic to calculate (default = mean)(options: mean, min, max)',  default='mean', nargs='?')
    parser.add_argument('-p', '--plev', help='flag to select a pressure level (defaults to None if not provided', nargs='?')

    
    args = parser.parse_args()


    if args.variable != None and args.frequency != None and args.timeslice != None and args.diri != None and args.diro != None:
        var=args.variable
        freq=args.frequency
        timeslice=args.timeslice
        diri=args.diri
        diro=args.diro
        area=args.area
        seas=args.seas
        stat=args.stat
        plev=args.plev

        #print("Start Rnet calc for year "+str(args.year)+" in "+str(args.diri)+" ...")
        #calc_rnet_tonetcdf(diri,freq,year)
        print("Start calc for "+str(args.timeslice)+" "+str(args.area)+" "+str(args.seas)+" "+str(args.stat)+" "+str(args.variable)+" in "+str(args.diro)+" plev="+str(args.plev)+" ...")

        #calc_ens_statistics_and_to_netcdf(var, freq, timeslice, diri, diro, area='global', seas='ANN', statistic='mean')
        calc_ens_statistics_and_to_netcdf(var, freq, timeslice,plev, diri, diro, area, seas, stat)      
        

if __name__ == '__main__':
    sys.exit(main())


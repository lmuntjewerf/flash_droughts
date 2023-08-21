import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

##############
## All variables/settings in the study

# Data
datadir = '/nobackup/users/wiel/rhine_myd/'
reanalyses = ['ERA5','JRA55','ERA20C','20CRv3']
models = ['CANESM2','CESM','CSIRO','ECEarth','GFDLCM3','GFDLESM2M','MPI']
rean_models = reanalyses.copy()
rean_models.extend(models)

# Multi-year drought settings/choices
spei_ref = [1961,2020]
CMD_spei_period = 6 # 6-SPEI
CMD_threshold = -1
CMD_min_months = 2
CMD_consec_years = 2
CMD_season = [6,11]  # months, inclusive
LDHD_spei_period = 12 # SPEI-12
LDHD_threshold = -1
LDHD_min_length = 12

# Figure settings 
model_ripf = {'CANESM2':'i1p1','CESM':'i1p1','CSIRO':'i1p1','ECEarth':'i1p1','GFDLCM3':'i1p1','GFDLESM2M':'i1p1','MPI':'i1850p3'}
model_names = {'ERA5':'ERA5','JRA55':'JRA-55','ERA20C':'ERA-20C','20CRv3':'20CRv3','CANESM2':'CanESM2','CESM':'CESM1-CAM5','CSIRO':'CSIRO-Mk3-6-0','ECEarth':'EC-Earth 2.3','GFDLCM3':'GFDL-CM3','GFDLESM2M':'GFDL-ESM2M','MPI':'MPI-ESM-LR'}
model_colors = {'ERA5':'black','JRA55':'black','ERA20C':'black','20CRv3':'black','CANESM2':'#a6cee3','CESM':'#1f78b4','CSIRO':'#33a02c','ECEarth':'#b2df8a','GFDLCM3':'#e31a1c','GFDLESM2M':'#fb9a99','MPI':'#fdbf6f'}
#model_colors = {'ERA5':'black','JRA55':'black','ERA20C':'black','20CRv3':'black','CANESM2':plt.cm.Set1(0),'CESM':plt.cm.Set1(1),'CSIRO':plt.cm.Set1(2),'ECEarth':plt.cm.Set1(3),'GFDLCM3':plt.cm.Set1(4),'GFDLESM2M':'gold','MPI':plt.cm.Set1(6)}
model_linestyles = {'ERA5':'-','JRA55':'--','ERA20C':'-.','20CRv3':':','CANESM2':'-','CESM':'-','CSIRO':'-','ECEarth':'-','GFDLCM3':'-','GFDLESM2M':'-','MPI':'-'}
start_years = {'ERA5':1950,'JRA55':1958,'ERA20C':1900,'20CRv3':1836,'CANESM2':1950,'CESM':1920,'CSIRO':1850,'ECEarth':1860,'GFDLCM3':1920,'GFDLESM2M':1950,'MPI':1850}
end_years = {'ERA5':2020,'JRA55':2020,'ERA20C':2010,'20CRv3':2015,'CANESM2':2100,'CESM':2100,'CSIRO':2100,'ECEarth':2100,'GFDLCM3':2100,'GFDLESM2M':2100,'MPI':2099}
ensmembers = {'ERA5':1,'JRA55':1,'ERA20C':1,'20CRv3':1,'CANESM2':50,'CESM':35,'CSIRO':30,'ECEarth':16,'GFDLCM3':20,'GFDLESM2M':30,'MPI':100}
alphabet = 'abcdefghijklmnopqrstuvwxyz'
var_longnames = {'tas':'temperature','pr':'precipitation','rsds':'incoming solar radiation','pet':'reference evapotranspiration','SPEI6':'SPEI-6','SPEI12':'SPEI-12','psl':'sea level pressure','wb':'precipitation deficit','swvl':'volumetric soil moisture','mrros':'Runoff','mrso':'soil moisture content'}
var_shortnames = {'tas':'TAS','pr':'PR','rsds':'RSDS','pet':'ET$_0$','SPEI6':'SPEI-6','SPEI12':'SPEI-12','psl':'PSL','wb':'ET$_0$ - PR','swvl':'SWVL','mrros':'MRROS','mrso':'MRSO'}
var_units = {'tas':'$^\circ$C','pr':'mm/month','rsds':'W/m$^2$','pet':'mm/month','SPEI6':'','SPEI12':'','psl':'Pa','wb':'mm/month','swvl':'10$^3$ m$^3$/m$^3$','mrros':'mm/month','mrso':'kg/m$^2$'}

##############
## Functions

def regions(name):
    if name == 'boxRhine':
        region = [2, 12, 46, 54]
        landonly = True
    elif name == 'global':
        region = [0, 360, -90, 90]
        landonly = False
    return region, landonly

def data_origensfile(model, var, ensmemb):
    """Open the (regridded) original ensemble data files, single member"""
    # open file
    if model in models:
        file = f"{datadir}{model}/rg/{var}_m_{model}_historical_rcp85_r{ensmemb}i*p*_*.nc"
    else:
        file = f"{datadir}{model}/rg/{var}_m_{model}.nc"
    print(file)
    ds = xr.open_mfdataset(file, decode_times=False)
    da = ds[var]    
    ds.close()
    # fix time
    new_time = pd.date_range(f"{start_years[model]}-01-01", periods=len(da.time), freq="M")
    da['time'] = new_time
    # fix units
    da = fix_units(da, var, model)
    # other random fixes
    if (model == 'CESM') & (ensmemb == 1):
        da = da.sel(time=slice('1920-01-1','2100-12-31'))
    # done
    return da

def data_fullensfiles(model, var):
    """Open the (regridded) original ensemble data files, all members"""
    # open file
    if model in models:
        file = f"{datadir}{model}/rg/{var}_m_{model}_historical_rcp85_r*i*p*_*.nc"
    else:
        file = f"{datadir}{model}/rg/{var}_m_{model}.nc"
    print(file)
    ds = xr.open_mfdataset(file, combine='nested', concat_dim='ensemble', decode_times=False)
    da = ds[var]    
    ds.close()
    # fix time
    new_time = pd.date_range(f"{start_years[model]}-01-01", periods=len(da.time), freq="M")
    da['time'] = new_time
    # fix units
    da = fix_units(da, var, model)
    # done
    return da

def data_tsensfile(model,var,region,file=None):
    """Open timeseries of ensemble data (output boxmean script)"""
    # open file
    if file is None:
        file = f"{datadir}{model}/{var}_m_{model}_{region}.nc"
    else:
        file = f"{datadir}{model}/{file}"
    print(file)
    ds = xr.open_mfdataset(file)
    da = ds[var]    
    ds.close()
    # done
    return da

def setlatlonweights(da, landonly, model):
    """Create weights for area-mean, possibly with land-sea contrast applied"""
    weights_lat = np.cos(da.lat/180*np.pi).values
    weights_latlon = np.reshape(np.repeat(weights_lat,len(da.lon)),(len(da.lat),len(da.lon)))
    if landonly: # Check paths for land-sea masks
        da_lsmask = xr.open_dataset(f"{datadir}{model}/rg/landsea_{model}.nc")['landsea']
        if da_lsmask.ndim == 3: 
            da_lsmask = da_lsmask[0]
        weights_latlon = weights_latlon * da_lsmask
    da_weights = xr.DataArray(weights_latlon,dims=['lat','lon'],coords=({'lat':da.lat,'lon':da.lon}),name='weights')
    return da_weights
def fldavg(da, dim=None, weights=None, ncfile=None):
    """
    Function to calculate weighted average over dimension. Mostly used for lat,lon averages. Can be used on xarray dataset through the apply function
    Example: ds.apply(fldavg,dim=['lat','lon'],weights=gridarea.cell_area)
    Courtesy: Folmer Krikken, KNMI (E: folmer.krikken@knmi.nl)
    """
    if weights is None and ncfile is None:
        return da.mean(dim)
    else:
        if da.notnull().any():
            total_weights = weights.where(da.notnull()).sum(dim=dim)
        else:
            total_weights = weights.sum(dim)
        return (da * weights).sum(dim) / total_weights

def fix_units(da,var,model): 
    """Set standard units"""  
    if var == 'pr':
        month_length = da.time.dt.days_in_month
        if model in models:
            # From mm/s to mm/month
            da_convert = da * 86400 * month_length
        elif model in ['ERA5','ERA20C']:
            # From m/d to mm/month
            da_convert = da * 1000 * month_length
        elif model == '20CRv3':
            # From kg/m2/s to mm/month
            da_convert = da * 86400 * month_length
        elif model == 'JRA55':
            # From mm/d to mm/month
            da_convert = da * month_length
        else:
            da_convert = da
        da_convert.attrs['units'] = 'mm/month'
        da_convert.name = var
    elif var == 'tas':
        if da.mean() > 200:
            # From K to C
            da_convert = da - 273.15
        else:
            da_convert = da
        da_convert.attrs['units'] = 'C'
    elif var == 'rsds':
        if model in ['ERA5','ERA20C']:
            # From J/s to W/m2
            da_convert = da / (24 * 60 * 60)
        else:
            da_convert = da
        da_convert.attrs['units'] = 'W m-2'
    elif var == 'mrros':
        month_length = da.time.dt.days_in_month
        if model in ['ERA5','ERA20C']:
            # From m/d to mm/month
            da_convert = da * 1000 * month_length
        else:
            da_convert = da
    elif var == 'swvl':
        if model in ['ERA5']:
            da_convert = da * 1000
        else:
            da_convert = da
    else:
        da_convert = da
    return da_convert

def Makkink_KNMI(T, Q):
    """Compute Makkink reference evapotranspiration (KNMI internal formula)
    Courtesy: Emma Daniels"""
    #1) verzadigde dampspanning tov water
    e_s = 6.107 * 10**(7.5 * (T / (237.3 + T)))   #[hPa]
    #2) verzadigde dampspanningsgradient tov water
    delta = ((7.5 * 237.3) /  (T + 237.3)**2) * np.log(10) * e_s   #[hPa/degC]
    #3) psychrometerconstante (afhankelijk van T)
    gamma = 0.646 + (0.0006 * T)   #[hPa/degC]
    #4) verdampingswarmte van water
    labda = 1000 * (2501 - (2.375 * T))   #[J/kg]
    #5) soortelijke massa van water
    rho = 1000   #[kg/m3]
    #verdamping
    month_length = T.time.dt.days_in_month
    Ev = ( (1000 * 0.65 * delta) / ((delta + gamma) * rho * labda) ) * (Q * (60 * 60 * 24 * month_length))   #[mm/month]
    # Karin: add attributes
    Ev.name = 'pet'
    Ev.attrs['long_name'] = 'Makkink reference evapotranspiration (KNMI)'
    Ev.attrs['standard_name'] = 'ref_evapotranspiration'
    Ev.attrs['units'] = 'mm/month'
    return Ev

def month_to_season(da,season):
    """Compute timeseries of seasonal means, note DJF is still JF...D"""
    if season in ['DJF','JJA','MAM','SON']:
        # temporal weights
        month_length = da.time.dt.days_in_month
        weights_mon = month_length.sel(time=da.time.dt.season==season)
        weights_season = weights_mon.groupby(weights_mon.time.dt.year).sum(dim='time')
        # mean (weighted by days in month)
        da_sel = da.sel(time=da.time.dt.season==season)
        da_season = (da_sel*weights_mon).groupby(da_sel.time.dt.year).sum(dim='time') / weights_season
    elif season == 'JJASON':
        # temporal weights
        month_length = da.time.dt.days_in_month
        weights_mon = month_length.where((da.time.dt.month>5) & (da.time.dt.month<12))
        weights_season = weights_mon.groupby(weights_mon.time.dt.year).sum(dim='time')
        # mean (weighted by days in month)
        da_sel = da.where((da.time.dt.month>5) & (da.time.dt.month<12))
        da_season = (da_sel*weights_mon).groupby(da_sel.time.dt.year).sum(dim='time') / weights_season
    else:
        # temporal weights
        month_length = da.time.dt.days_in_month
        weights_mon = month_length
        weights_season = weights_mon.groupby(weights_mon.time.dt.year).sum(dim='time')
        # mean (weighted by days in month)
        da_sel = da
        da_season = (da_sel*weights_mon).groupby(da_sel.time.dt.year).sum(dim='time') / weights_season
    # done
    return da_season

def CMD_identify_MeteoDroughts(da):
    """CMD identification step 1: identify all years with a meteorological drought"""
    # find all points under threshold
    count = da.where(da.values<CMD_threshold,99)
    count = count.where(count==99,1)
    count = count.where(count==1,0)
    count = xr.DataArray(count,dims=('time'),coords=({'time':da.time}))
    # number of points per year
    count_y = count.groupby(count.time.dt.year).sum(dim='time')
    # only years with minimum number/length
    years = count_y.where(count_y>=CMD_min_months).dropna(dim='year').year.values
    return years
def CMD_identify_ConsecutiveDroughts(years,ensmemb):
    """CMD identification step 2: find all consecutive years with a meteorological drought"""
    diff = years[1:]-years[:-1]
    couples = np.where(diff<2)[0]
    tmp = -10
    sel_no = 0
    consec_years = []
    for i in couples:
        if tmp != i-1:
            if sel_no != 0:
                consec_years.append(array)
            # NEW multi-year event!
            # store ensemble member
            array = [ensmemb]
            # story year 1 and year 2
            array.append(years[i])
            array.append(years[i+1])
            sel_no = sel_no + 1
            tmp = i
        else:
            # store year 3+++
            if years[i+1] != array[-1]:
                array.append(years[i+1])
            tmp = i
    if 'array' in locals():
        consec_years.append(array)
    return  consec_years
def CMD_writetofile(array,file_out): 
    """CMD identification step 3: write to txt file"""
    # aantal droogtes = aantal regels
    # lengte droogtes = lengte regel
    for i in range(len(array)):
        for j in array[i]:
            file_out.write(f"{int(j)} ")
        file_out.write('\n')
def CMD_readfile(filename):
    # read file
    file_in = open(filename)
    lines = file_in.readlines()
    file_in.close()
    # define arrays
    ensemble = []
    years = []
    length = []
    # fill arrays
    for line in lines:
        ensemble.append(int(line.split(' ')[0]))
        tmp = []
        for i in range(1,len(line[:-2].split())):
            tmp.append(int(line.split(' ')[i]))
        years.append(tmp)
        length.append(len(line[:-2].split())-1)
    # create dataframe
    data = {'ensemble':ensemble,'years':years,'length':length}
    df = pd.DataFrame(data)
    return df

def LDHD_identify(da,ensmemb,file_out):
    """LDHD identification: find them, write to file"""
    # find all points under threshold
    count = da.where(da>LDHD_threshold,99)
    count = count.where(count==99,0)
    count = count.where(count==0,1)
    # determine drought length
    tmp_len = 0
    tmp_yr = 0
    for i in range(np.where(count==0)[0][0],len(count)):
        if count[i] == 1: 
            tmp_len = tmp_len + 1
            if tmp_yr == 0:
                tmp_yr = da[i].time.dt.year.values
                tmp_mon = da[i].time.dt.month.values
        else:
            if tmp_len > LDHD_min_length-1:
                file_out.write(f"{ensmemb}\t{tmp_yr}-{str(tmp_mon).zfill(2)}\t{tmp_len}\n")
            tmp_len = 0
            tmp_yr = 0
    # check if there is a LDHD at the end of the run (add for count, not for duration)
    if tmp_len > LDHD_min_length-1:
        file_out.write(f"{ensmemb}\t{tmp_yr}-{str(tmp_mon).zfill(2)}\t9999\n")
def LDHD_readfile(filename):
    # read file
    file_in = open(filename)
    lines = file_in.readlines()
    file_in.close()
    # define arrays
    ensemble = []
    start = []
    year = []
    month = []
    length = []
    # fill arrays
    for i in range(len(lines)):
        line = lines[i][:-1].split('\t')
        ensemble.append(int(line[0]))
        start.append(line[1])
        year.append(int(line[1][0:4]))
        month.append(int(line[1][5:]))
        if int(line[2]) != 9999:
            length.append(int(line[2]))
        else:
            length.append(np.nan)
    # create dataframe
    data = {'ensemble':ensemble,'start':start,'year':year,'month':month,'length':length}
    df = pd.DataFrame(data)
    return df
    
def MYD_limitperiod(df,period,type='CMD'):
    period = np.arange(period[0],period[1]+1)
    if type == 'CMD':
        ensemble = []
        years = []
        length = []
        for i_e in range(len(df)):
            if df.iloc[i_e].years[0] in period:
                ensemble.append(df.iloc[i_e].ensemble)
                years.append(df.iloc[i_e].years)
                length.append(df.iloc[i_e].length)
        data = {'ensemble':ensemble,'years':years,'length':length}
    elif type == 'LDHD':
        ensemble = []
        start = []
        year = []
        month = []
        length = []
        for i_e in range(len(df)):
            if df.iloc[i_e].year in period:
                ensemble.append(df.iloc[i_e].ensemble)
                start.append(df.iloc[i_e].start)
                year.append(df.iloc[i_e].year)
                month.append(df.iloc[i_e].month)
                length.append(df.iloc[i_e].length)
        data = {'ensemble':ensemble,'start':start,'year':year,'month':month,'length':length}
    df_period = pd.DataFrame(data)
    return df_period

def detrend_poly(da,model,base_period,order,type):
    n_ens = ensmembers[model]
    da_month = da.where(da.time.dt.month==1).dropna(dim='time')
    yrs = da_month.time.dt.year.values
    yrs_1D = np.tile(yrs,n_ens)
    da_detrend = np.full(da.shape,np.nan)
    # detrend by month
    if type == 'bymonth':
        for month in np.arange(1,13):
            # an array of only this month
            da_month = da.where(da.time.dt.month==month).dropna(dim='time')
            # make 1D
            da_month_1D = np.reshape(da_month.values,n_ens*len(yrs))
            # fit polynomial
            trend_function = np.poly1d(np.polyfit(yrs_1D,da_month_1D,order))
            da_month_trend = trend_function(yrs)
            da_month_trend = xr.DataArray(da_month_trend,dims=['time'],coords=({'time':da_month.time}))
            # detrend
            da_month_detrend = da_month - da_month_trend + da_month_trend.sel(time=slice(f"{base_period[0]}-01-01",f"{base_period[1]}-12-31")).mean(dim='time')
            # save in data-array
            da_detrend[:,month-1::12] = da_month_detrend
            """# possible: check
            plt.plot(yrs,da_month[2],color='grey',linewidth=.5)
            plt.plot(yrs,da_month_trend,color='red')
            plt.plot(yrs,da_month_detrend[2],color='black',linewidth=.5)
            plt.show()
            import sys; sys.exit()"""
    # detrend by year
    elif type == 'byyear':
        # annual mean
        da_yr = da.groupby(da.time.dt.year).mean(dim='time')
        # make 1D
        da_yr_1D = np.reshape(da_yr.values,n_ens*len(yrs))
        yrs_1D = np.tile(da_yr.year,n_ens)
        # fit polynomial
        trend_function = np.poly1d(np.polyfit(yrs_1D,da_yr_1D,order))
        da_yr_trend = trend_function(yrs)
        da_yr_trend = xr.DataArray(da_yr_trend,dims=['year'],coords=({'year':da_yr.year}))
        # detrend
        da_detrend = da.groupby(da.time.dt.year) - da_yr_trend + da_yr_trend.sel(year=slice(f"{base_period[0]}",f"{base_period[1]}")).mean(dim='year')
        # maybe: check check double check
        """plt.plot(yrs,da_yr[0],color='black',linewidth=.5)
        plt.plot(yrs,da_yr_trend,color='red')
        plt.show()
        yrs_m = np.arange(yrs[0],yrs[-1]+1,1/12)
        plt.plot(yrs_m,da[2],color='black',linewidth=.5)
        plt.plot(yrs_m,da_detrend[2],color='red',linewidth=.5)
        plt.show()"""
    # xarray data array
    da_detrend = xr.DataArray(da_detrend,dims=['ensemble','time'],coords=[da.ensemble,da.time],name=da.name)
    da_detrend.attrs = da.attrs
    if type == 'bymonth':
        da_detrend.attrs['detrended'] = f"{order}-order polynomial, base period {base_period[0]}-{base_period[1]}, detrended by month"
    elif type == 'byyear':
        da_detrend.attrs['detrended'] = f"{order}-order polynomial, base period {base_period[0]}-{base_period[1]}, detrended by year"
    return da_detrend

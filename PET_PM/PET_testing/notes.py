
from /perm/nklm/mambaforge/lib/python3.10/site-packages/pyet/./utils.py 

def get_index(df):
    """Method to return the index of the input data.

    """
    try:
        index = pandas.DatetimeIndex(df.index)
    except AttributeError:
        index = pandas.DatetimeIndex(df.time)
    return index

from /perm/nklm/mambaforge/lib/python3.10/site-packages/pyet/rad_utils.py

def calc_rad_net(tmean, rn=None, rs=None, lat=None, n=None, nn=None, tmax=None,
                 tmin=None, rhmax=None, rhmin=None, rh=None, elevation=None,
                 rso=None, a=1.35, b=-0.35, ea=None, albedo=0.23, as1=0.25,
                 bs1=0.5, kab=None):
    """Net radiation [MJ m-2 d-1].

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    kab: float, optional
        coefficient derived from as1, bs1 for estimating clear-sky radiation
        [degrees].

    Returns
    -------
    float/pandas.Series/xarray.DataArray, optional containing the calculated
        net shortwave radiation
        Notes
    -----
    Based on equation 40 in :cite:t:`allen_crop_1998`.

    """
    import numpy as np

    if rn is not None:
        rn = check_rad(rn)
        return rn
    else:
        if rs is None:
            rs = calc_rad_sol_in(n, lat, as1=as1, bs1=bs1, nn=nn)
        rns = calc_rad_short(rs=rs, lat=lat, n=n, nn=nn, albedo=albedo,
                             as1=as1, bs1=bs1)  # [MJ/m2/d]
        print(np.shape(rs))
        print(np.shape(tmean))
        print(np.shape(tmax))
        print(np.shape(tmin))
        print(np.shape(rhmax))
        print(np.shape(rhmin))
        print(np.shape(rh))
        print(np.shape(elevation))
        print(np.shape(lat))
        print(np.shape(rso))
        print(np.shape(a))
        print(np.shape(b))
        print(np.shape(ea))
        print(np.shape(kab))
        rnl = calc_rad_long(rs=rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, rso=rso, a=a, b=b,
                            ea=ea, kab=kab)  # [MJ/m2/d]
        #print(np.shape(rns))
        #print(np.shape(rnl))
        #print(len(np.shape(rnl)))
        if len(np.shape(rnl))==3:
            rnl['lat']=rns['lat']
        rn = rns- rnl
        print(np.shape(rn))
        #rn = check_rad(rn)
        #print(np.shape(rn))
        return rn

# ---> verder met calc_rad_long, daar zit het probleem. 
    if rso is None:
        #tindex = get_index(rs)
        tindex = get_index(tmean)


# het probleem is niet met de 23:00 uur stamp.. want zowel voor de 1D resample.mean() en voor de cdo daymean is het een probleem. 
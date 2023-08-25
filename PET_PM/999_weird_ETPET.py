
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


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
    #Ev = ( (1000 * 0.65 * delta) / ((delta + gamma) * rho * labda) ) * (Q * (60 * 60 * 24 * month_length))   #[mm/month]
    Ev = ( (1000 * 0.65 * delta) / ((delta + gamma) * rho * labda) ) * (Q * (60 * 60 * 24 ))   #[mm/day?]
    # Karin: add attributes
    Ev.name = 'pet'
    Ev.attrs['long_name'] = 'Makkink reference evapotranspiration (KNMI)'
    Ev.attrs['standard_name'] = 'ref_evapotranspiration'
    Ev.attrs['units'] = 'mm/month'
    return Ev







datadir = '/perm/nkkw/P2_flashdroughts/link_project/meteodata_ERA5/'


"""
da_tas = xr.open_dataset(f"{datadir}tas_Rhine.nc").tas
da_rsds = xr.open_dataset(f"{datadir}rsds_Rhine.nc").rsds
da_makk = Makkink_KNMI(da_tas,da_rsds)
da_makk.name = 'makk'
ds_out = da_makk.to_dataset()
ds_out.to_netcdf(f"{datadir}makk_Rhine.nc")

da_tas = xr.open_dataset(f"{datadir}tas_Rhine_detrend.nc").tas
da_rsds = xr.open_dataset(f"{datadir}rsds_Rhine_detrend.nc").rsds
da_makk = Makkink_KNMI(da_tas,da_rsds)
da_makk.name = 'makk'
da_makk.attrs['detrend'] = 'computed using detrended TAS and RSDS'
ds_out = da_makk.to_dataset()
ds_out.to_netcdf(f"{datadir}makk_Rhine_detrend.nc")

import pdb; pdb.set_trace()
"""



da_et = xr.open_dataset(f"{datadir}et_Rhine.nc").et
da_pet = xr.open_dataset(f"{datadir}pet_Rhine.nc").pet
da_makk = xr.open_dataset(f"{datadir}makk_Rhine.nc").makk


year = 2020
da_et = da_et.sel(time=slice(f"{year}-01-01",f"{year}-12-31"))
da_pet = da_pet.sel(time=slice(f"{year}-01-01",f"{year}-12-31"))
da_makk = da_makk.sel(time=slice(f"{year}-01-01",f"{year}-12-31"))

da_esi = da_et/da_pet
da_esi_makk = da_et / da_makk

print()
print("\t\tERA5\t\t\t\tMakkink")
print("time\t\tET\tPET\tESI\t\tPET\tESI")
for i in [0,1,2,3,345,346,347]:
  print(f"{da_et.time[i].dt.day.values}-{da_et.time[i].dt.month.values}\t\t{np.round(da_et[i].values,4)}\t{np.round(da_pet[i].values,4)}\t {np.round(da_esi[i].values,2)}\t\t{np.round(da_makk[i].values,4)}\t{np.round(da_esi_makk[i].values,2)}")
print()

fig = plt.figure(figsize=(9,4))
plt.plot(da_et,color='firebrick',label='ERA5_ET')
plt.plot(da_pet,color='deepskyblue',label='ERA5_PET')
plt.plot(da_makk,color='royalblue',label='Makkink_PET')
plt.legend()

fig = plt.figure(figsize=(9,4))
plt.plot(da_esi,color='deepskyblue',label='ERA5_ESI')
plt.plot(da_esi_makk,color='royalblue',label='Makkink_ESI')
plt.ylim(-.5,5)
plt.legend()

plt.show()




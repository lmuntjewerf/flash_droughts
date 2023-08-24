

FD_index = 'SMI-7'
FD_window = 28
FD_jump = 2.0
#FD_startthreshold = 0
FD_endthreshold = -1.5

basin = 'Rhine'

plot_indices = True; plot_ind = ['SPI-28','SPEI-28','ESI-14','SMI-7']
plot_meteo = True; plot_met = ['pr','et','mrsos','tas']

#-----------------------------------------

import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import re
import matplotlib.pyplot as plt

datadir = '/perm/nkkw/P2_flashdroughts/link_project/'
# figure_out_dir = 'link_project/figures/'
figure_out_dir = '/perm/nklm/Px_flashdroughts/figures'
fontsize=10
alphabet = 'abcdefghij'
maanden = ['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
longnames = {'pr':'Precipitation','tas':'2m temperature','rsds':'Solar radiation','et':'Actual and potential evapotranspiration','mrsos':'Surface soil moisture'}
units = {'pr':'mm/day','tas':'$^\circ$C','rsds':'W/m2','et':'mm/day','mrsos':'m3/m3','pet':'mm/day'}

def func_plot_indices(variables=['SPI21','SPEI21','ESI14','SMI7']):
    # make figure
    plt.figure(figsize=(7,7))
    plt.subplots_adjust(left=0.10, right=0.90, top=0.97, bottom=0.04, wspace=.5, hspace=.5)
    if len(variables) > 4:
        print('ERROR: too many indices to plot')
        import sys; sys.exit()
    for i,var in enumerate(variables):
        ax = plt.subplot(4,1,i+1)
        # open data
        ds = xr.open_dataset(f"{datadir}composites/{var}_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
        da_comp = ds[var]
        ds.close()
        # plot
        xaxis = da_comp.time.values
        plt.plot(xaxis,da_comp.mean(dim='event'),color='steelblue',linewidth=2,zorder=10)
        #for j in range(10):
        #    plt.plot(xaxis,da_comp[j],color='grey',linewidth=.75,alpha=.6,zorder=5)
        #plt.fill_between(xaxis,da_comp.quantile(0,dim='event'),da_comp.quantile(1,dim='event'),color='steelblue',alpha=.3,linewidth=0)
        #plt.fill_between(xaxis,da_comp.quantile(0.25,dim='event'),da_comp.quantile(.75,dim='event'),color='steelblue',alpha=.4,linewidth=0)
        # annotation
        if 'FD_index' in globals():
            if var == FD_index:
                plt.plot([0,0],[-4,4],color='firebrick',linewidth=1,linestyle='-')
                if ('FD_startthreshold' in globals()) and ('FD_window' in globals()):
                    plt.plot([-5,5],np.repeat(FD_startthreshold,2),color='firebrick',linewidth=2,zorder=7)
                if ('FD_endthreshold' in globals()) and ('FD_window' in globals()):
                    plt.plot([FD_window-5,FD_window+5],np.repeat(FD_endthreshold,2),color='firebrick',linewidth=2,zorder=7)
                if ('FD_jump' in globals()) and ('FD_window' in globals()):
                    plt.plot([1.5,2.5,2,2,1.5,2.5],[-.8,-.8,-.8,-.8-FD_jump,-.8-FD_jump,-.8-FD_jump],color='firebrick',linewidth=1,zorder=7)
        if 'FD_window' in globals():
            plt.fill_between([0,FD_window],[-4-4],[4,4],color='firebrick',alpha=.2,zorder=0,linewidth=0)
        plt.title(f"{alphabet[i]}) {var}",fontsize=fontsize,loc='left')
        if i == 0:
            plt.title(f"FD criteria: {FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}",fontsize=fontsize,loc='right')
        plt.ylabel(f"{var} [{da_comp.units}]",fontsize=fontsize-1,color='steelblue')
        plt.gca().xaxis.grid(True)
        plt.plot([-70,100],[0,0],color='black',linewidth=0.75,zorder=0)
        plt.ylim(-3,3)
        plt.yticks(np.arange(-3,3.1,1),fontsize=fontsize-1)
        plt.xticks([-60,-45,-30,-15,0,15,30,45,60,75,90],fontsize=fontsize-1)
        plt.xlim(-60,60)

def func_plot_meteo(variables=['pr','et','mrsos','tas'],plot_met_anom=True):
    # make figure
    plt.figure(figsize=(7,7))
    plt.subplots_adjust(left=0.10, right=0.90, top=0.97, bottom=0.04, wspace=.5, hspace=.5)
    if len(variables) > 4:
        print('ERROR: too many variables to plot')
        import sys; sys.exit()
    for i,var in enumerate(variables):
        ax = plt.subplot(4,1,i+1)
        # open data
        if not plot_met_anom:
            ds = xr.open_dataset(f"{datadir}composites/{var}_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            typename = ''
        else:
            ds = xr.open_dataset(f"{datadir}composites/{var}_anom_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            typename = 'anomalies'
        da_comp = ds[var]
        ds.close()
        if var == 'et':
            if not plot_met_anom:
                ds_extra = xr.open_dataset(f"{datadir}composites/pet_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            else:
                ds_extra = xr.open_dataset(f"{datadir}composites/pet_anom_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            da_extra = ds_extra['pet']
        elif var == 'tas':
            if not plot_met_anom:
                ds_extra = xr.open_dataset(f"{datadir}composites/rsds_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            else:
                ds_extra = xr.open_dataset(f"{datadir}composites/rsds_anom_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.nc")
            da_extra = ds_extra['rsds']
        # plot
        xaxis = da_comp.time.values
        plt.plot(xaxis,da_comp.mean(dim='event'),color='steelblue',linewidth=2,zorder=10)
        #for j in range(10):
        #    plt.plot(xaxis,da_comp[j],color='grey',linewidth=.75,alpha=.6,zorder=5)
        #plt.fill_between(xaxis,da_comp.quantile(0,dim='event'),da_comp.quantile(1,dim='event'),color='steelblue',alpha=.3,linewidth=0)
        #plt.fill_between(xaxis,da_comp.quantile(0.25,dim='event'),da_comp.quantile(.75,dim='event'),color='steelblue',alpha=.4,linewidth=0)
        if 'FD_window' in globals():
            ylim = ax.get_ylim()
            plt.fill_between([0,FD_window],np.repeat(np.min(da_comp.values)-10,2),np.repeat(np.max(da_comp.values)+10,2),color='firebrick',alpha=.2,zorder=0,linewidth=0)
        if var != 'tas':
            plt.title(f"{alphabet[i]}) {longnames[var]} {typename}",fontsize=fontsize,loc='left')
        else:
            plt.title(f"{alphabet[i]}) {longnames[var]} & {longnames['rsds']} {typename}",fontsize=fontsize,loc='left')
        if i == 0:
            plt.title(f"FD criteria: {FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}",fontsize=fontsize,loc='right')
        plt.ylabel(f"{var.upper()} [{units[var]}]",color='steelblue',fontsize=fontsize-1)
        plt.gca().xaxis.grid(True)
        plt.plot([-70,100],[0,0],color='black',linewidth=0.75,zorder=0)
        plt.ylim(ylim[0],ylim[1])
        plt.yticks(fontsize=fontsize-1)
        plt.xticks([-60,-45,-30,-15,0,15,30,45,60,75,90],fontsize=fontsize-1)
        if var in ['et','tas']:
            ax2=ax.twinx()
            ax2.plot(xaxis,da_extra.mean(dim='event'),color='deepskyblue',linewidth=2)
            if var == 'et':
                plt.ylabel(f"PET [{units['pet']}]",color='deepskyblue',fontsize=fontsize-1)
                ylim2 = ax2.get_ylim()
                ax.set_ylim(np.min([ylim[0],ylim2[0]]),np.max([ylim[1],ylim2[1]]))
                ax2.set_ylim(np.min([ylim[0],ylim2[0]]),np.max([ylim[1],ylim2[1]]))
            elif var == 'tas':
                plt.ylabel(f"RSDS [{units['rsds']}]",color='deepskyblue',fontsize=fontsize-1)
                plt.plot([-70,100],[0,0],color='black',linewidth=0.75,linestyle='--',zorder=0)
            plt.yticks(fontsize=fontsize-1)
        plt.xlim(-30,45)

def running_mean(da,rm):
    da_rm = np.convolve(da.values, np.ones(rm)/rm, mode='same')
    da_rm = xr.DataArray(da_rm,coords=da.coords,name=da.name)
    da_rm.attrs = da.attrs
    return da_rm

if plot_indices:
    func_plot_indices(plot_ind)
    plt.savefig(f"{figure_out_dir}/comp_ind_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.png")
if plot_meteo:
    func_plot_meteo(plot_met,True)
    plt.savefig(f"{figure_out_dir}/comp_metA_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.png")
    func_plot_meteo(plot_met,False)
    plt.savefig(f"{figure_out_dir}/comp_met_{basin}_{FD_index}_w{FD_window}_j{int(FD_jump*10)}_e{int(abs(FD_endthreshold)*10)}.png")
#plt.show()




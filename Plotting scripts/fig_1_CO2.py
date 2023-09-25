import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import xesmf as xe
import regionmask
from cmip6_preprocessing.regionmask import merged_mask

#%%
year_start=2015
year_end1=2100
year_end2=2100
RUN1='control' # List of runs: control, hosing_05, 
RUN2='hosing_05'
quality=300 # quality when figure is saved in dpi
a=1         # amount of years averaged over
tm=a
save_fig='yes'

A = ([2025, 2040, 2055, 2070, 2085, 2100])
LW=4
FS=18

VARS='TMCO2','-'
UNIT= '[ppm]','-'
CONV= 1,12.01*1e-8,1e-3,1e6*28.966/44
DESCR='CO$_2$','-'

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
area_gr=load_var1['areacello'][:,:].compute().squeeze()

lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Prepare regridder
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(area_gn, ds_out, 'bilinear',periodic=True)       
    
#%% Location data
data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_126/'
data_snel2='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_126/'
data_snel3='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_585/'
data_snel4='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_585/' 

n=0
#for n in range (len(VARS)):
var=VARS[n]
unit=UNIT[n]
conv=CONV[n]
descr=DESCR[n]
  
#%% Call on datasets (might need to change name of dataset)
load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()
    
load_var1 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR2=load_var1[var][:,:,:].compute().squeeze()

plt.plot(VAR2['time'])

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR4=load_var1[var].compute().squeeze()

#%% Regrid data
VAR1_gr = VAR1
VAR2_gr = VAR2
VAR3_gr = VAR3
VAR4_gr = VAR4

lat=VAR1_gr.lat
lon=VAR1_gr.lon

RE = 6.371e6  # [m] Earth radius

dy = 2*np.pi*RE*(lat[1]-lat[0]).values/360                              # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360    # grid size in x-direction

VAR_lon=(VAR1_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var1=VAR_lat/5.14e18*1.0e6 * 28.966 / 44.0

VAR_lon2=(VAR2_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat2=(VAR_lon2[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var2=VAR_lat2/5.14e18*1.0e6 * 28.966 / 44.0

VAR_lon=(VAR3_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var3=VAR_lat/5.14e18*1.0e6 * 28.966 / 44.0

VAR_lon2=(VAR4_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat2=(VAR_lon2[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var4=VAR_lat2/5.14e18*1.0e6 * 28.966 / 44.0
    
plt_VAR1 = pd.DataFrame(plt_var1)
rolling_windows = plt_VAR1.rolling(tm)
PLT_VAR1=np.squeeze(rolling_windows.mean())

plt_VAR2 = pd.DataFrame(plt_var2)
rolling_windows = plt_VAR2.rolling(tm)
PLT_VAR2=np.squeeze(rolling_windows.mean())
    
plt_VAR1 = pd.DataFrame(plt_var3)
rolling_windows = plt_VAR1.rolling(tm)
PLT_VAR3=np.squeeze(rolling_windows.mean())

plt_VAR2 = pd.DataFrame(plt_var4)
rolling_windows = plt_VAR2.rolling(tm)
PLT_VAR4=np.squeeze(rolling_windows.mean())

#%%
t1=np.arange(2016,2100.5,1)  
t2=t1

datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'
 
fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,PLT_VAR3,linewidth=LW,color='tab:blue',linestyle='solid',label='CTL-585')
ax.plot(t2,PLT_VAR4,linewidth=LW,color='tab:orange',linestyle='solid',label='HOS-585')
ax.plot(t1,PLT_VAR1,linewidth=LW,color='tab:blue',linestyle='dashed',label='CTL-126')
ax.plot(t1,PLT_VAR2[:],linewidth=LW,color='tab:orange',linestyle='dashed',label='HOS-126')

plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel(descr+' ' + unit,fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)
plt.xlim([2020,2100])

if save_fig == 'yes':
    plt.savefig(datadir+'fig_1_c.png', format='png', dpi=quality,bbox_inches='tight')

fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,PLT_VAR4-PLT_VAR3,linewidth=LW,color='tab:green',label='585')
ax.plot(t1,np.array(PLT_VAR2)[:]-PLT_VAR1,linewidth=LW,color='tab:green',linestyle='dashed',label='126')
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('$\Delta$'+descr+' ' + unit,fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.xlim([2020,2100])
plt.legend(fontsize=FS-4)
plt.ylim([-4,6])

if save_fig == 'yes':
    plt.savefig(datadir+'fig_1_f.png', format='png', dpi=quality,bbox_inches='tight')

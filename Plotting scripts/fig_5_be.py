import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import xesmf as xe
import cmocean as cm
import regionmask
from cmip6_preprocessing.regionmask import merged_mask

#%%
year_start=2015
year_end=2100
year_end1=2100
RUN1='control'
RUN2='hosing_05' # List of runs: control, hosing_05, 
quality=300 # quality when figure is saved in dpi
a=5         # amount of years averaged over
save_fig='yes'
years=np.arange(year_start,year_end+1,1)
A = ([2025, 2040, 2055, 2070, 2085, 2100])

VARS='FG_CO2', '-'   
UNIT1='[kg C s$^{-1}$]', '-'
UNIT2='[kg C]', '[-]'
UNIT3='[kg C s$^{-1}$ m$^{-1}$]', '-'


DESCR='Gas exchange','-'
scen='126' # 126 for 5b or 585 for 5e

if scen == '126':
    sub = 'b'
    style = 'dashed'
elif scen == '585':
    sub = 'e'
    style = 'solid'

#%%
tm = 5
def mov_mean(data1, data2, data3):
    Lv11 = pd.DataFrame(data1)
    rolling_windows = Lv11.rolling(tm)
    var1=np.squeeze(rolling_windows.mean())
    
    Lv11 = pd.DataFrame(data2)
    rolling_windows = Lv11.rolling(tm)
    var2=np.squeeze(rolling_windows.mean())
    
    Lv11 = pd.DataFrame(data3)
    rolling_windows = Lv11.rolling(tm)
    var3=np.squeeze(rolling_windows.mean())
    
    return var1,var2,var3

#%%
LW = 4
FS = 20

#%%
def plot_2(data1,data2,data3,data4,data5,var,unit,delta):
    plt.plot(t,np.cumsum(data1*365*86400*1e-12),linewidth = LW,linestyle = style, label = 'Atlantic')
    plt.plot(t,np.cumsum(data2*365*86400*1e-12),linewidth = LW,linestyle = style, label = 'Indian')
    plt.plot(t,np.cumsum(data3*365*86400*1e-12),linewidth = LW,linestyle = style, label = 'Pacific')
    plt.plot(t,np.cumsum(data4*365*86400*1e-12),linewidth = LW,linestyle = style, label = 'Southern')
    plt.plot(t,np.cumsum(data5*365*86400*1e-12),linewidth = LW,linestyle = style, label = 'Arctic')
    
    plt.xlabel('Time',fontsize=FS-2)
    if delta == 1:
        text = '$\Delta$' + var + ' ' + unit
    else:
        text = var + ' ' + unit
    plt.ylabel(text,fontsize = FS-2)
    plt.xticks(A,fontsize=FS-4)
    plt.yticks(fontsize=FS-4)
    plt.xlim([2020,2100])
    plt.grid()
    plt.legend(fontsize=FS-6)
    
#%%
data1='/Users/daan/Documents/Snellius_stuff/Repository/Data/'   # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
area_gr=load_var1['areacello'][:,:].compute().squeeze()

#%% Define basins
basins = regionmask.defined_regions.natural_earth.ocean_basins_50
mask = merged_mask(basins,area_gr.squeeze())

lon = area_gr['lon']
lat = area_gr['lat']
mask = merged_mask(basins,area_gr.squeeze())

#%% Define masks
# Atlantic
masked1 = np.logical_or(np.logical_or(mask==0,mask==0),mask==1)
masked1=np.array(masked1)*np.ones((np.shape(mask)))
masked1[:55,:]=0
masked1[157:,:]=0
masked1[masked1==0]='nan'
masked1=np.array(masked1)

# Indian
masked2 = np.logical_or(np.logical_or(mask==5,mask==5),mask==5)
masked2=np.array(masked2)*np.ones((np.shape(mask)))
masked2[:55,:]=0
masked2[masked2==0]='nan'
masked2=np.array(masked2)

# Pacific
masked3 = np.logical_or(np.logical_or(mask==2,mask==3),mask==2)
masked3=np.array(masked3)*np.ones((np.shape(mask)))
masked3[:55,:]=0
masked3[masked3==0]='nan'
Masked3=np.array(masked3)

#%% Prepare regridder
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(area_gn, ds_out, 'bilinear',periodic=True)       
    
#%% Gas exchange plots (Figure 1a, b and c)
if scen == '585':
    data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_585/'   # Location of dataset(s) 
    data_snel2='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_585/'   # Location of dataset(s) 
else: 
    data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_126/'   # Location of dataset(s) 
    data_snel2='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_126/'   #
    
n=0
var=VARS[n]
unit1=UNIT1[n]
unit2=UNIT2[n]
unit3=UNIT3[n]
unit4=unit3
descr=DESCR[n]

print(var)
   
#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end1)+'_0'+str(scen)+'.nc' )
if load_var2[var][:,0,0].size == 192 or load_var2[var][:,0,0].size == 384:
    VAR3=load_var2[var][:,:,:].compute().squeeze()
    VAR2=VAR3.transpose('time',...)
else:
    VAR2=load_var2[var][:,:,:].compute().squeeze()

print(np.shape(VAR2))
load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end1)+'_0'+str(scen)+'.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

#%%
t = np.arange(2016,2100.5,1)
t = t[:]

#%%
VAR1_gr = regridder(VAR1)
VAR_ctl=np.roll(VAR1_gr,-180)   

VAR2_gr = regridder(VAR2)
VAR_hos=np.roll(VAR2_gr,-180) 

#%%
masked4=masked3
masked4[:55,:]=VAR_ctl[0,:55,:]/VAR_ctl[0,:55,:]
masked4[55:,:]='nan'
masked4=np.array(masked4)

#%%
masked5=masked3
masked5[158:,:]=VAR_ctl[0,158:,:]/VAR_ctl[0,158:,:]
masked5[:158,:]='nan'

#%%
VAR_ctl_1 = VAR_ctl*masked1
VAR_ctl_2 = VAR_ctl*masked2
VAR_ctl_3 = VAR_ctl*Masked3
VAR_ctl_4 = VAR_ctl*masked4
VAR_ctl_5 = VAR_ctl*masked5

VAR_hos_1 = VAR_hos*masked1
VAR_hos_2 = VAR_hos*masked2
VAR_hos_3 = VAR_hos*Masked3
VAR_hos_4 = VAR_hos*masked4
VAR_hos_5 = VAR_hos*masked5

#%%
VAR2_ctl_1 = np.nansum(np.nansum(VAR_ctl_1*np.array(area_gr),axis=1),axis=1)
VAR2_ctl_2 = np.nansum(np.nansum(VAR_ctl_2*np.array(area_gr),axis=1),axis=1)
VAR2_ctl_3 = np.nansum(np.nansum(VAR_ctl_3*np.array(area_gr),axis=1),axis=1)
VAR2_ctl_4 = np.nansum(np.nansum(VAR_ctl_4*np.array(area_gr),axis=1),axis=1)
VAR2_ctl_5 = np.nansum(np.nansum(VAR_ctl_5*np.array(area_gr),axis=1),axis=1)

VAR2_hos_1 = np.nansum(np.nansum(VAR_hos_1*np.array(area_gr),axis=1),axis=1)
VAR2_hos_2 = np.nansum(np.nansum(VAR_hos_2*np.array(area_gr),axis=1),axis=1)
VAR2_hos_3 = np.nansum(np.nansum(VAR_hos_3*np.array(area_gr),axis=1),axis=1)
VAR2_hos_4 = np.nansum(np.nansum(VAR_hos_4*np.array(area_gr),axis=1),axis=1)
VAR2_hos_5 = np.nansum(np.nansum(VAR_hos_5*np.array(area_gr),axis=1),axis=1)

dVAR2_1 = VAR2_hos_1-VAR2_ctl_1
dVAR2_2 = VAR2_hos_2-VAR2_ctl_2
dVAR2_3 = VAR2_hos_3-VAR2_ctl_3
dVAR2_4 = VAR2_hos_4-VAR2_ctl_4
dVAR2_5 = VAR2_hos_5-VAR2_ctl_5

#%%
[var2_ctl_1,var2_hos_1,dvar2_1] = mov_mean(VAR2_ctl_1,VAR2_hos_1,dVAR2_1)
[var2_ctl_2,var2_hos_2,dvar2_2] = mov_mean(VAR2_ctl_2,VAR2_hos_2,dVAR2_2)
[var2_ctl_3,var2_hos_3,dvar2_3] = mov_mean(VAR2_ctl_3,VAR2_hos_3,dVAR2_3)
[var2_ctl_4,var2_hos_4,dvar2_4] = mov_mean(VAR2_ctl_4,VAR2_hos_4,dVAR2_4)
[var2_ctl_5,var2_hos_5,dvar2_5] = mov_mean(VAR2_ctl_5,VAR2_hos_5,dVAR2_5)

#%%
datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'

#%%
fig = plt.figure(figsize=(7, 5))
plot_2(dvar2_1,dvar2_2,dvar2_3,dvar2_4,dvar2_5,descr,'[PgC]',1)
plt.ylim([-10,5.5])
if save_fig == 'yes':
    plt.savefig(datadir+'fig_5'+str(sub)+'.png', format='png', dpi=quality,bbox_inches='tight')

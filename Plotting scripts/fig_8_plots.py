import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import xesmf as xe
import cmocean as cm
import gsw

#%%
year_start=2015
year_end_126=2100
year_end_585=2100

RUN1='control'
RUN2='hosing_05' # List of runs: control, hosing_05, 
quality=300 # quality when figure is saved in dpi
a=5         # amount of years averaged over
save_fig='yes'
years=np.arange(year_start+1,year_end_585+1,1)

VARS= 'NBP','FG_CO2'
UNIT = '[kg C m$^{-2}$]','[kg C m$^{-2}$]'
DESCR = 'NBP','Gas exchange'

cmap1 = 'RdYlGn'
cmap2 = 'BrBG_r'

scen = '585'
lat1 = -90
lat2 = 90
lon1 = -179-60
lon2 = 179-60

reg = 'glob'

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
area_gr=load_var1['areacello'][:,:].compute().squeeze()

#%% Prepare regridder
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(area_gn, ds_out, 'bilinear',periodic=True)       
    
#%% Location data
data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_126/'
data_snel2='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_126/'
data_snel3='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_585/'
data_snel4='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_585/' 

#%%
def subplot(data1,data2,data3,data4,data5,data6,scen):
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(0))
    #ax.set_extent([lon1+60, lon2+60, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='110m',zorder=6)
    
    im1=plt.pcolormesh(lon1,lat1,data2-data1,vmin=-4,vmax=4,transform=ccrs.PlateCarree(),cmap=cmap1)
    im2=plt.pcolormesh(lon2,lat2,data4-data3,vmin=-1.25,vmax=1.25,transform=ccrs.PlateCarree(),cmap=cmap2)
    plt.contour(lon3,lat3,-data5,levels = [-0.25],colors='purple',linewidths=(4.5))
    plt.contour(lon3,lat3,data6,levels = [0.25],colors='purple',linewidths=(4.5))
    ax.set_title( str(scen),fontsize=FS+5)
    
    cax1 = fig.add_axes([0.15, 0.21, 0.325, 0.02])
    cbar1 = fig.colorbar(im1, cax=cax1, orientation = 'horizontal')
    
    # Colorbar specifics
    cbar1.ax.set_xlabel('$\Delta$Land exchange [kg C m$^{-2}$]', fontsize=16) 
    cbar1.ax.set_yticks(fontsize=16)
    cbar1.ax.tick_params(labelsize=16)
    cbar1.ax.xaxis.offsetText.set_fontsize(16)
    
    cax2 = fig.add_axes([0.55, 0.21, 0.325, 0.02])
    cbar2 = fig.colorbar(im2, cax=cax2, orientation = 'horizontal')
    cbar2.ax.set_xlabel('$\Delta$Ocean exchange [kg C m$^{-2}$]', fontsize=16) 
    cbar2.ax.set_yticks(fontsize=16)
    cbar2.ax.tick_params(labelsize=16)
    cbar2.ax.xaxis.offsetText.set_fontsize(16)
     
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,zorder=5)
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}

#%%
NN = 0
var=VARS[NN]


#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

#%% Regrid data
a = VAR1[0,:,:].where(np.isnan(VAR1[0,:,:]),1)

VAR1_gr1 = np.sum(VAR1,axis=0)*86400*365*a
VAR2_gr1 = np.sum(VAR2,axis=0)*86400*365*a
VAR3_gr1 = np.sum(VAR3,axis=0)*86400*365*a
VAR4_gr1 = np.sum(VAR4,axis=0)*86400*365*a

#%%
lat1 = VAR1.lat
lon1 = VAR1.lon
    
#%%
NN = 1
var=VARS[NN]

#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

#%%
VAR1_gr = regridder(VAR1)
VAR1_gr=np.roll(VAR1_gr,-180)*86400*365 

VAR2_gr = regridder(VAR2)
VAR2_gr=np.roll(VAR2_gr,-180)*86400*365

VAR3_gr = regridder(VAR3)
VAR3_gr=np.roll(VAR3_gr,-180)*86400*365

VAR4_gr = regridder(VAR4)
VAR4_gr=np.roll(VAR4_gr,-180)*86400*365 

#%%
b = np.where(np.isnan(VAR1_gr[0,:,:]),VAR1_gr[0,:,:],1)

VAR1_gr2 = np.sum(VAR1_gr,axis=0)*b
VAR2_gr2 = np.sum(VAR2_gr,axis=0)*b
VAR3_gr2 = np.sum(VAR3_gr,axis=0)*b
VAR4_gr2 = np.sum(VAR4_gr,axis=0)*b

lat2 = area_gr.lat
lon2 = area_gr.lon

#%%
var='ICEFRAC'

#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

#%%
b = np.where(np.isnan(VAR1[0,:,:]),VAR1[0,:,:],1)

VAR1_gr3 = np.mean(VAR1[25:35,:,:],axis=0)*b
VAR2_gr3 = np.mean(VAR2[25:35,:,:],axis=0)*b
VAR3_gr3 = np.mean(VAR3[25:35,:,:],axis=0)*b
VAR4_gr3 = np.mean(VAR4[25:35,:,:],axis=0)*b

lat3 = VAR1.lat
lon3 = VAR1.lon

#%%
var='TREFHT'

#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

#%%
VAR1_gr = np.nanmean(VAR1[-5:,:,:],axis=0)
VAR2_gr = np.nanmean(VAR2[-5:,:,:],axis=0)
VAR3_gr = np.nanmean(VAR3[-5:,:,:],axis=0)
VAR4_gr = np.nanmean(VAR4[-5:,:,:],axis=0)

VAR1_gr4 = np.nanmean(VAR1_gr[:,:],axis=1)
VAR2_gr4 = np.nanmean(VAR2_gr[:,:],axis=1)
VAR3_gr4 = np.nanmean(VAR3_gr[:,:],axis=1)
VAR4_gr4 = np.nanmean(VAR4_gr[:,:],axis=1)

lat4 = VAR1.lat
lon4 = VAR1.lon
     
#%% Plotting variables
FS=30

datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'

scen = 'SSP 5-8.5'

fig = plt.figure(figsize=(20, 16))
subplot(VAR3_gr1,VAR4_gr1,VAR3_gr2,VAR4_gr2,VAR3_gr3,VAR4_gr3,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10b.png', format='png', dpi=quality,bbox_inches='tight')


#%%
scen = 'SSP 1-2.6'

fig = plt.figure(figsize=(20, 16))
subplot(VAR1_gr1,VAR2_gr1,VAR1_gr2,VAR2_gr2,VAR1_gr3,VAR2_gr3,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10a.png', format='png', dpi=quality,bbox_inches='tight')
    

#%%
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap

fig = plt.figure(figsize=(10, 0.5))

ax = fig.add_axes([0, 0, 1, 1])
#ax.set_axis_off()

plt.xticks([])
plt.yticks([])
# create a collection with a rectangle for each year
col = PatchCollection([
    Rectangle((y, 0), 1, 1)
    for y in range(-90, 90 + 1)
])

# set data, colormap and color limits

col.set_array(VAR2_gr4-VAR1_gr4)
col.set_cmap('RdBu_r')
col.set_clim([-2.5,2.5])
ax.add_collection(col)
ax.set_ylim(0, 1)
ax.set_xlim(-90, 90 + 1)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10a_temp.png', format='png', dpi=quality,bbox_inches='tight')

#%%
fig = plt.figure(figsize=(10, 0.5))

ax = fig.add_axes([0, 0, 1, 1])
#ax.set_axis_off()

# create a collection with a rectangle for each year
col = PatchCollection([
    Rectangle((y, 0), 1, 1)
    for y in range(-90, 90 + 1)
])
plt.xticks([])
plt.yticks([])

# set data, colormap and color limits

col.set_array(VAR4_gr4-VAR3_gr4)
col.set_cmap('RdBu_r')
col.set_clim([-2.5,2.5])
ax.add_collection(col)
ax.set_ylim(0, 1)
ax.set_xlim(-90, 90 + 1)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10b_temp.png', format='png', dpi=quality,bbox_inches='tight')



#%%
d1_5 = VAR4_gr1-VAR3_gr1
d1_5a = d1_5.where(d1_5>0)
d1_5b = d1_5.where(d1_5<0)

d1_1 = VAR2_gr1-VAR1_gr1
d1_1a = d1_1.where(d1_5>0)
d1_1b = d1_1.where(d1_5<0)

d2_5 = VAR4_gr2-VAR3_gr2
d2_5a = np.where(d2_5>0,d2_5,np.nan)
d2_5b = np.where(d2_5<0,d2_5,np.nan)

d2_1 = VAR2_gr2-VAR1_gr2
d2_1a = np.where(d2_5>0,d2_1,np.nan)
d2_1b = np.where(d2_5<0,d2_1,np.nan)

Mask1a = np.isnan(d1_5a)
Mask2a = np.isnan(d2_5a)
Mask1b = np.isnan(d1_5b)
Mask2b = np.isnan(d2_5b)
    

#%%
scen = 'Difference (negative anomalies)'

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(0))
#ax.set_extent([lon1+60, lon2+60, lat1, lat2], ccrs.PlateCarree())
    
ax.coastlines(resolution='110m',zorder=6)

im1=plt.pcolormesh(lon1,lat1,d1_5b-d1_1b,vmin=-4,vmax=4,transform=ccrs.PlateCarree(),cmap='PiYG')
im2=plt.pcolormesh(lon2,lat2,d2_5b-d2_1b,vmin=-1.25,vmax=1.25,transform=ccrs.PlateCarree(),cmap='PuOr_r')
ax.set_title( str(scen),fontsize=22)

plt.contour(lon1,lat1,Mask1b,levels = [0],colors='black',linewidths=(0.5))
plt.contour(lon2,lat2,Mask2b,levels = [np.nan],colors='black',linewidths=(0.5))
 
# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,zorder=5)
gl.xlabel_style = {'size': 16}
gl.ylabel_style = {'size': 16}

cax1 = fig.add_axes([0.15, 0.17, 0.7, 0.035])
cbar1 = fig.colorbar(im1, cax=cax1, orientation = 'horizontal')

# Colorbar specifics
cbar1.ax.set_xlabel('$\Delta$Land exchange [kg C m$^{-2}$]', fontsize=16) 
cbar1.ax.set_yticks(fontsize=16)
cbar1.ax.tick_params(labelsize=16)
cbar1.ax.xaxis.offsetText.set_fontsize(16)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10c.png', format='png', dpi=quality,bbox_inches='tight')

#%%
scen = 'Difference (positive anomalies)'

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(0))
#ax.set_extent([lon1+60, lon2+60, lat1, lat2], ccrs.PlateCarree())
    
ax.coastlines(resolution='110m',zorder=6)

im1=plt.pcolormesh(lon1,lat1,d1_5a-d1_1a,vmin=-4,vmax=4,transform=ccrs.PlateCarree(),cmap='PiYG')
im2=plt.pcolormesh(lon2,lat2,d2_5a-d2_1a,vmin=-1.25,vmax=1.25,transform=ccrs.PlateCarree(),cmap='PuOr_r')
ax.set_title( str(scen),fontsize=22)

plt.contour(lon1,lat1,Mask1a,levels = [0],colors='black',linewidths=(0.5))
plt.contour(lon2,lat2,Mask2a,levels = [np.nan],colors='black',linewidths=(0.5))
 
# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,zorder=5)
gl.xlabel_style = {'size': 16}
gl.ylabel_style = {'size': 16}
    
cax2 = fig.add_axes([0.15, 0.17, 0.7, 0.035])
cbar2 = fig.colorbar(im2, cax=cax2, orientation = 'horizontal')
cbar2.ax.set_xlabel('$\Delta$Ocean exchange [kg C m$^{-2}$]', fontsize=16) 
cbar2.ax.set_yticks(fontsize=16)
cbar2.ax.tick_params(labelsize=16)
cbar2.ax.xaxis.offsetText.set_fontsize(16)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_10d.png', format='png', dpi=quality,bbox_inches='tight')

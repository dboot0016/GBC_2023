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
quality=300 # quality when figure is saved in dpi
save_fig='yes'

var= 'NBP'

cmap = 'Greys'

scen = '585'
lat1 = 30
lat2 = 90
lon1 = -90
lon2 = 30

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
    
#%% Gas exchange plots (Figure 1a, b and c)
data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/'   # Location of dataset(s) 
FS = 18

#%%
def subplot(data1):
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
    ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND,zorder=4)
    ax.coastlines(resolution='110m',zorder=6)

    im=plt.pcolormesh(lon,lat,data1,transform=ccrs.PlateCarree(),vmin=0,vmax=0.9,cmap=cmap)
    ax.set_title('Hosing region',fontsize=FS)
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    ax.text(-4000000,4100000,'40$^{\circ}$N',fontsize=14)
    ax.text(-4000000,5200000,'50$^{\circ}$N',fontsize=14)
    ax.text(-4000000,6200000,'60$^{\circ}$N',fontsize=14)
    ax.text(-4000000,7100000,'70$^{\circ}$N',fontsize=14)
    ax.text(-4000000,8000000,'80$^{\circ}$N',fontsize=14)
    
    ax.text(-2800000,2700000,'270$^{\circ}$E',fontsize=14)
    ax.text(-500000,2700000,'300$^{\circ}$E',fontsize=14)
    ax.text(2000000,2700000,'330$^{\circ}$E',fontsize=14)
    ax.text(5000000,2700000,'0$^{\circ}$E',fontsize=14)
    ax.text(7200000,2700000,'30$^{\circ}$E',fontsize=14)


#%% Call on datasets (might need to change name of dataset)
load_var1 = xr.open_dataset(f'{data_snel1}/FW_BSSP585_BPRPcmip6_hosing_final.nc' )
VAR1=load_var1['FW_liquid'][0,:,:].compute().squeeze()

#%%
VAR1_gr = regridder(VAR1)
VAR1_gr=np.roll(VAR1_gr,-180)

lat = area_gr.lat
lon = area_gr.lon

#%%
datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'

fig = plt.figure(figsize=(7, 5))
subplot(VAR1_gr/VAR1_gr)
if save_fig == 'yes':
    plt.savefig(datadir+'fig_S1.png', format='png', dpi=quality,bbox_inches='tight')

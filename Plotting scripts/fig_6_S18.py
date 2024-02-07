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

VARS= 'NBP','COL_FIRE_CLOSS'
UNIT = '[kg C m$^{-2}$]','[kg C m$^{-2}$]'
DESCR = 'NBP','Fire loss','CH$_4$ oxidation'
CMAP = (['RdYlGn','YlOrBr'])
FIG_NR = '8','S17'

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
def subplot(data1,data2,i,scen):
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(0))
    ax.set_extent([lon1+60, lon2+60, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='50m',zorder=6)
    ax.add_feature(cfeature.OCEAN,zorder=4)
    
    if i == 0:
        im=plt.pcolormesh(lon,lat,np.nansum(data1[:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Uptake (CTL-' + str(scen)+')',fontsize=FS)
    
    elif i == 1:
        im=plt.pcolormesh(lon,lat,np.nansum(data2[:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Uptake (HOS-' + str(scen)+')',fontsize=FS)
        
    elif i == 2:
        im=plt.pcolormesh(lon,lat,np.nansum(data2[:,:,:],axis=0)-np.nansum(data1[:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference (HOS-CTL: ' + str(scen)+')',fontsize=FS)

    # Colorbar specifics
    cbar=plt.colorbar(im,orientation='horizontal',) 
    if i == 0 or i == 1:
        cbar.ax.set_xlabel(descr+' ' + unit, fontsize=16)
    elif i == 2:
        cbar.ax.set_xlabel('$\Delta$'+descr+' ' + unit, fontsize=16)
        
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,zorder=5)

    if lon1 == (-179-60) and not var == 'ICEFRAC':
        ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
        ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
        ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
        ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
        ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

        ax.text(-10500000,-10000000,'180$^{\circ}$E',fontsize=14)
        ax.text(-6000000,-10000000,'270$^{\circ}$E',fontsize=14)
        ax.text(-1000000,-10000000,'0$^{\circ}$E',fontsize=14)
        ax.text(3000000,-10000000,'90$^{\circ}$E',fontsize=14)
        ax.text(7000000,-10000000,'180$^{\circ}$E',fontsize=14)

#%%
for NN in range(len(VARS)):
    var=VARS[NN]
    unit=UNIT[NN]
    descr=DESCR[NN]
    cmap=CMAP[NN]
    fig_nr=FIG_NR[NN]
    
    print(var)
    
          
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
    VAR1_gr = VAR1*86400*365
    VAR2_gr = VAR2*86400*365
    VAR3_gr = VAR3*86400*365
    VAR4_gr = VAR4*86400*365
    
    lat = VAR1.lat
    lon = VAR1.lon
         
    #%% Plotting variables
    FS=20
        
    #%%
    datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'
    
    VMN1=([-10,0])
    VMX1=([10,10])
    
    vmn1 = VMN1[NN]
    vmx1 = VMX1[NN]
    
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR3_gr,VAR4_gr,0,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR3_gr,VAR4_gr,1,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'e.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    VMN1=([-10,0])
    VMX1=([10,10])
    
    vmn1 = VMN1[NN]
    vmx1 = VMX1[NN]
    
    scen = '126'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR1_gr,VAR2_gr,0,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR1_gr,VAR2_gr,1,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'b.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    VMN1=([-6,-4])
    VMX1=([6,4])
    
    vmn1 = VMN1[NN]
    vmx1 = VMX1[NN]
    
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR3_gr,VAR4_gr,2,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'f.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR1_gr,VAR2_gr,2,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'c.png', format='png', dpi=quality,bbox_inches='tight')

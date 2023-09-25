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
year_start1=2016
year_start2=2096
year_end1=2020
year_end2=2100

RUN1='control'
RUN2='hosing_05' # List of runs: control, hosing_05, 
quality=300 # quality when figure is saved in dpi
a=5         # amount of years averaged over
save_fig='yes'

VARS='WVEL','DIC_flux','PO4_flux' 
UNIT = '[m day$^{-1}$]','[mol m$^{-2}$ day$^{-1}$]','[mol m$^{-2}$ day$^{-1}$]'
DESCR = 'w$_{(150m)}$','DIC upwelling (150m)','PO$_4$ upwelling (150m)'
CMAP = (['PuOr','PuOr','PuOr'])
FIG_NR = 'S10','S14','S15'

NN = 0

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
def subplot(data1,data2,i,exp,scen):
    if var == 'ICEFRAC':
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    else:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
        ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='50m',zorder=5)
    #if var == 'NBP':
     #   ax.add_feature(cfeature.OCEAN)
        
    if not descr == 'SAT'  and not var =='TAUX' and not var =='TAUY' and not 'NBP' and not 'PRECT':
        ax.add_feature(cfeature.LAND,zorder=4)
        
    ax.add_feature(cfeature.LAND,zorder=4)
    
    if i == 0:
        im=plt.pcolormesh(lon,lat,data1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Average 2016-2020 ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
        
    elif i == 1:
        im=plt.pcolormesh(lon,lat,data2,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Average 2096-2100 ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
        
    # Colorbar specifics
    cbar=plt.colorbar(im,orientation='horizontal',) 
    cbar.ax.set_xlabel(descr+' ' + unit, fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    
    if lon1 == (-179-60) and not var == 'ICEFRAC':
        ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
        ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
        ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
        ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
        ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)
        
        ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
        ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
        ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
        ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

#%%
for NN in range(len(VARS)):
    var=VARS[NN]
    unit=UNIT[NN]
    descr=DESCR[NN]
    cmap=CMAP[NN]
    fig_nr=FIG_NR[NN]
    
    print(var)
       
    #%% Call on datasets (might need to change name of dataset)
    if var == 'WVEL':
        load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR21=load_var2[var].compute().squeeze()*1e-2*86400
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR11=load_var1[var].compute().squeeze()*1e-2*86400
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR41=load_var2[var].compute().squeeze()*1e-2*86400
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR31=load_var1[var].compute().squeeze()*1e-2*86400
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR22=load_var2[var].compute().squeeze()*1e-2*86400
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR12=load_var1[var].compute().squeeze()*1e-2*86400
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR42=load_var2[var].compute().squeeze()*1e-2*86400
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR32=load_var1[var].compute().squeeze()*1e-2*86400
        
    elif var == 'DIC_flux' or var == 'PO4_flux':
        vara = 'WVEL'
        if var == 'DIC_flux':
            varb = 'DIC1'
            varc = 'DIC1'
        elif var == 'PO4_flux':
            varb = 'PO4'
            varc = 'PO41'
            
        load_var2 = xr.open_dataset(f'{data_snel2}/'+vara+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR2a=load_var2[vara].compute().squeeze()*1e-2*86400
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+vara+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR1a=load_var1[vara].compute().squeeze()*1e-2*86400
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+vara+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR4a=load_var2[vara].compute().squeeze()*1e-2*86400
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+vara+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR3a=load_var1[vara].compute().squeeze()*1e-2*86400
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+varc+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR2b=load_var2[varb].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+varc+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR1b=load_var1[varb].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+varc+'_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR4b=load_var2[varb].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+varc+'_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR3b=load_var1[varb].compute().squeeze()
        
        VAR11 = VAR1a*VAR1b
        VAR21 = VAR2a*VAR2b
        VAR31 = VAR3a*VAR3b
        VAR41 = VAR4a*VAR4b
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+vara+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR2a=load_var2[vara].compute().squeeze()*1e-2*86400
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+vara+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR1a=load_var1[vara].compute().squeeze()*1e-2*86400
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+vara+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR4a=load_var2[vara].compute().squeeze()*1e-2*86400
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+vara+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR3a=load_var1[vara].compute().squeeze()*1e-2*86400
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+varc+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR2b=load_var2[varb].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+varc+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR1b=load_var1[varb].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+varc+'_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR4b=load_var2[varb].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+varc+'_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR3b=load_var1[varb].compute().squeeze()
        
        VAR12 = VAR1a*VAR1b
        VAR22 = VAR2a*VAR2b
        VAR32 = VAR3a*VAR3b
        VAR42 = VAR4a*VAR4b
    
    #%% Regrid data
    VAR1_gr1 = regridder(VAR11.mean('time'))
    VAR1_gr1=np.roll(VAR1_gr1,-180)   
    
    VAR2_gr1 = regridder(VAR21.mean('time'))
    VAR2_gr1=np.roll(VAR2_gr1,-180) 
    
    VAR3_gr1 = regridder(VAR31.mean('time'))
    VAR3_gr1=np.roll(VAR3_gr1,-180)   
    
    VAR4_gr1 = regridder(VAR41.mean('time'))
    VAR4_gr1=np.roll(VAR4_gr1,-180) 
    
    VAR1_gr2 = regridder(VAR12.mean('time'))
    VAR1_gr2=np.roll(VAR1_gr2,-180)   
    
    VAR2_gr2 = regridder(VAR22.mean('time'))
    VAR2_gr2=np.roll(VAR2_gr2,-180) 
    
    VAR3_gr2 = regridder(VAR32.mean('time'))
    VAR3_gr2=np.roll(VAR3_gr2,-180)   
    
    VAR4_gr2 = regridder(VAR42.mean('time'))
    VAR4_gr2=np.roll(VAR4_gr2,-180) 
    
    lat = area_gr.lat
    lon = area_gr.lon
         
    #%% Plotting variables
    FS=20
        
    #%%
    datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'
    
    Vmn1 = ([-0.45,-0.3,-6e-4])
    Vmx1 = ([0.45,0.3,6e-4])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    scen = '585'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR3_gr1,VAR3_gr2,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    scen = '585'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR4_gr1,VAR4_gr2,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'_extra_HOS_585_2016.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR3_gr1,VAR3_gr2,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'_extra_CTL_585_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR4_gr1,VAR4_gr2,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'_extra_HOS_585_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    Vmn1 = ([-0.45,-0.3,-6e-4])
    Vmx1 = ([0.45,0.3,6e-4])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    scen = '126'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR1_gr1,VAR1_gr2,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR2_gr1,VAR2_gr2,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'extra_HOS_126_2016.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    scen = '126'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR1_gr1,VAR1_gr2,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'extra_CTL_126_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR2_gr1,VAR2_gr2,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_'+str(fig_nr)+'extra_HOS_126_2096.png', format='png', dpi=quality,bbox_inches='tight')

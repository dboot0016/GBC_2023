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

VARS= 'ICEFRAC','FG_CO22','FG_CO2','TEMP','SALT','photoC_TOT','pH_3D','HMXL','DIC','ALK'
UNIT = '[-]','[kg C m$^{-2}$]','[g C m$^{-2}$ yr$^{-1}$]','[$^{\circ}$C]','[g/kg]','[mol C m$^{-2}$ s$^{-1}$]','[-]','[m]','[mol m$^{-2}$]','[mol m$^{-2}$]'
DESCR = 'Ice fraction','Gas exchange','Gas exchange','SST','SSS','NPP','pH','MLD','DIC (0-150m)','ALK (0-150m)'
CMAP = ([cm.cm.thermal,'BrBG_r','seismic','seismic',cm.cm.ice,'RdYlGn','BrBG_r',cm.cm.thermal,cm.cm.haline,cm.cm.ice,cm.cm.algae_r,cm.cm.matter_r,'RdYlBu','RdYlBu','RdYlBu','RdYlBu','RdYlBu','RdYlBu','PiYG','viridis',cm.cm.deep,'RdYlBu_r','RdYlBu_r','RdYlBu_r'])
CONV = ([1,1,1e3*86400*365,1,1,1,1,1,1,1])

NN = -1

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
    
#%% Gas exchange plots (Figure 1a, b and c)
data_snel1='/Users/daan/Documents/Snellius_stuff/control_126/'   # Location of dataset(s) 
data_snel2='/Users/daan/Documents/Snellius_stuff/hosing_126/'   # Location of dataset(s) 

data_snel3='/Users/daan/Documents/Snellius_stuff/control/'   # Location of dataset(s) 
data_snel4='/Users/daan/Documents/Snellius_stuff/hosing/'

#%%
def subplot(data1,a,scen):
    if scen == 'SSP1-2.6':
        ice1 = ice_126_c
        ice2 = ice_126_h
        
    elif scen == 'SSP5-8.5':
        ice1 = ice_585_c
        ice2 = ice_585_h
        
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NearsidePerspective(-40,55,9600000))
    #ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='110m',zorder=5)
        
    ax.add_feature(cfeature.LAND,zorder=4)  
    if var == 'FG_CO22':
        im=plt.pcolormesh(lon,lat,(np.sum(data1[:,:,:],axis=0))*86400*365,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        plt.contour(lon_ice,lat_ice,-np.mean(ice1[-5:,:,:],axis=0),[-0.15],colors='black',transform=ccrs.PlateCarree())
        plt.contour(lon_ice,lat_ice,np.mean(ice2[-5:,:,:],axis=0),[0.15],colors='black',transform=ccrs.PlateCarree())    
        ax.set_title(str(scen),fontsize=FS)
    else:
        if a == 1:
            im=plt.pcolormesh(lon,lat,np.mean(data1[-5:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
            plt.contour(lon_ice,lat_ice,-np.mean(ice1[-5:,:,:],axis=0),[-0.15],colors='black',transform=ccrs.PlateCarree())
            plt.contour(lon_ice,lat_ice,np.mean(ice2[-5:,:,:],axis=0),[0.15],colors='black',transform=ccrs.PlateCarree())    
            ax.set_title(str(scen)+': 2096-2100',fontsize=FS)
            
        else:
            im=plt.pcolormesh(lon,lat,np.mean(data1[20:25,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
            plt.contour(lon_ice,lat_ice,-np.mean(ice1[20:25,:,:],axis=0),[-0.15],colors='black',transform=ccrs.PlateCarree())
            plt.contour(lon_ice,lat_ice,np.mean(ice2[20:25,:,:],axis=0),[0.15],colors='black',transform=ccrs.PlateCarree())    
            ax.set_title(str(scen)+': 2036-2040',fontsize=FS)

    # Colorbar specifics
    cbar=plt.colorbar(im,orientation='horizontal',pad=0.06) 
    cbar.ax.set_xlabel('$\Delta$'+descr+' ' + unit, fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)

#%%
def dens(salt,temp,lev):
    CT1=gsw.CT_from_pt(salt,temp)
    p1=gsw.p_from_z(lev,salt['lat'])
    
    rho=gsw.density.rho(salt,CT1,p1)

    return rho

#%%
var = VARS[0]
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
ice_126_h=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
ice_126_c=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
ice_585_h=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
ice_585_c=load_var1[var][:,:,:].compute().squeeze()

lon_ice = ice_126_h.lon
lat_ice = ice_126_h.lat

#%%
for NN in range(len(VARS)-1):#len(VARS)-6):
    var=VARS[NN+1]
    unit=UNIT[NN+1]
    descr=DESCR[NN+1]
    cmap=CMAP[NN+1]
    conv=CONV[NN+1]
    
    print(var)
    
          
    #%% Call on datasets (might need to change name of dataset)
    if var == 'FG_CO22':
        varn = 'FG_CO2'
        load_var2 = xr.open_dataset(f'{data_snel2}/'+varn+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
        VAR2=load_var2[varn][:,:,:].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+varn+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
        VAR1=load_var1[varn][:,:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+varn+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
        VAR4=load_var2[varn][:,:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+varn+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
        VAR3=load_var1[varn][:,:,:].compute().squeeze()
        
    else:
        load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
        VAR2=load_var2[var][:,:,:].compute().squeeze()*conv
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
        VAR1=load_var1[var][:,:,:].compute().squeeze()*conv
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
        VAR4=load_var2[var][:,:,:].compute().squeeze()*conv
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
        VAR3=load_var1[var][:,:,:].compute().squeeze()*conv
    
    var_126 = VAR2 - VAR1
    var_585 = VAR4 - VAR3

#%% Regrid data
    VAR_126 = regridder(var_126)
    VAR_126=np.roll(VAR_126,-180)   
    
    VAR_585 = regridder(var_585)
    VAR_585=np.roll(VAR_585,-180) 
     
    lat = area_gr.lat
    lon = area_gr.lon
         
    #%% Plotting variables
    FS=20
        
    #%%
    datadir = '/Users/daan/Documents/Snellius_stuff/Article 5 fig/'
    
    Vmn1=([-1.6,-14,-4.1,-5.1,-1.6e-7,-0.08,-28,-40,-50])
    
    vmn1=Vmn1[NN]
    vmx1=-vmn1
    
    scen = 'SSP5-8.5'
    a = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(5, 6))
    subplot(VAR_585,a,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_draft5_na_arc_'+var+'_contour_2096_'+str(585)+'.png', format='png', dpi=quality)
    
    scen = 'SSP1-2.6'
    fig = plt.figure(figsize=(5, 6))
    subplot(VAR_126,a,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_draft5_na_arc_'+var+'_contour_2096_'+str(126)+'.png', format='png', dpi=quality)
    
    scen = 'SSP5-8.5'
    a = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(5, 6))
    subplot(VAR_585,a,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_draft5_na_arc_'+var+'_contour_2036_'+str(585)+'.png', format='png', dpi=quality)
    
    scen = 'SSP1-2.6'
    fig = plt.figure(figsize=(5, 6))
    subplot(VAR_126,a,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_draft5_na_arc_'+var+'_contour_2036_'+str(126)+'.png', format='png', dpi=quality)
    
  

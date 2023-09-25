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

VARS= 'TREFHT','PRECT','TAUX','ICEFRAC_N','ICEFRAC_S','NBP','FG_CO2','TEMP','SALT','Strat','photoC_TOT','pH_3D','DIC','ALK','PO4'
UNIT = '[$^{\circ}$C]','[mm day$^{-1}$]','[N m$^{-2}$]','[-]','[-]','[kg C m$^{-2}$ s$^{-1}$]','[kg C m$^{-2}$ s$^{-1}$]','[$^{\circ}$C]','[g/kg]','[kg m$^{-3}$]','[mol C m$^{-2}$ s$^{-1}$]','[-]','[mol m$^{-2}$]','[mol m$^{-2}$]','[mol m$^{-2}$]'
DESCR = 'SAT','Precipitation',r'$\tau_{x}$','Ice fraction','Ice fraction','NBP','Gas exchange','SST','SSS','Stratification','NPP','pH','DIC (0-150m)','ALK (0-150m)','PO$_4$ (0-150m)'
CMAP = ([cm.cm.thermal,'BrBG_r','seismic',cm.cm.ice,cm.cm.ice,'RdYlGn','BrBG_r',cm.cm.thermal,cm.cm.haline,cm.cm.ice,cm.cm.algae_r,'viridis','RdYlBu_r','RdYlBu_r',cm.cm.algae])
FIG_NR = '2','S2','S3','S4','S5','9','4','S7','S8','S9','S6','S16','S12','S13','S11'

lat1 = -90
lat2 = 90
lon1 = -179-60
lon2 = 179-60

#%% Load in areas for regridder
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

#%% Plotting function
def subplot(data1,i,exp,scen):
    if var == 'ICEFRAC' and reg == 'north':
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    elif var == 'ICEFRAC' and reg == 'south':
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
    else:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
        ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='50m',zorder=5)
    if var == 'NBP':
        ax.add_feature(cfeature.OCEAN,zorder=4)
        
    if not descr == 'SAT'  and not var =='TAUX' and not var =='TAUY' and not 'NBP' and not 'PRECT':
        ax.add_feature(cfeature.LAND,zorder=4)
        
    ax.add_feature(cfeature.LAND,zorder=4)
    
    if i == 0:
        im=plt.pcolormesh(lon,lat,data1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Average 2016-2020 ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
        
    elif i == 1:
        im=plt.pcolormesh(lon,lat,data1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cmap)
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

#%% Function to compute density
def dens(salt,temp,lev):
    CT1=gsw.CT_from_pt(salt,temp)
    p1=gsw.p_from_z(lev,salt['lat'])
    
    rho=gsw.density.rho(salt,CT1,p1)

    return rho

#%% Select variable, units, description, colormap, and figure number (for saving)
for NN in range(len(VARS)):
    
    if VARS[NN] == 'ICEFRAC_N':
        var='ICEFRAC'
        cmap=CMAP[NN]
        reg = 'north'
    
    elif VARS[NN] == 'ICEFRAC_S':
        var='ICEFRAC'
        reg = 'south'
    
    else:
        var=VARS[NN]
        reg = 'global'
    
    unit=UNIT[NN]
    descr=DESCR[NN]
    cmap=CMAP[NN]
    fig_nr=FIG_NR[NN]

    print(var)
       
    #%% Call on datasets
    if var == 'PRECT':
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECC_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARa1=load_var2['PRECC'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'PRECC_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARb1=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECL_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARa2=load_var2['PRECL'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/PRECL_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARb2=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR21 = (VARa1+VARa2)*1e3*86400
        VAR11 = (VARb1+VARb2)*1e3*86400
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECC_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARa1=load_var2['PRECC'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'PRECC_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARb1=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECL_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARa2=load_var2['PRECL'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/PRECL_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARb2=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR41 = (VARa1+VARa2)*1e3*86400
        VAR31 = (VARb1+VARb2)*1e3*86400
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECC_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARa12=load_var2['PRECC'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'PRECC_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARb12=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECL_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARa22=load_var2['PRECL'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/PRECL_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARb22=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR22 = (VARa12+VARa22)*1e3*86400
        VAR12 = (VARb12+VARb22)*1e3*86400
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECC_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARa12=load_var2['PRECC'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'PRECC_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARb12=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECL_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARa22=load_var2['PRECL'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/PRECL_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARb22=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR42 = (VARa12+VARa22)*1e3*86400
        VAR32 = (VARb12+VARb22)*1e3*86400
        
    elif var == 'Strat':
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR21 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR11 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR41 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR31 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR22 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR12 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR42 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR32 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
    else:
        load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR21=load_var2[var][:,:].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
        VAR11=load_var1[var][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR41=load_var2[var][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
        VAR31=load_var1[var][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR22=load_var2[var][:,:].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
        VAR12=load_var1[var][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR42=load_var2[var][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
        VAR32=load_var1[var][:,:].compute().squeeze()

    #%% Regrid data if necessary
    
    if var == 'TREFHT':
        VAR11_gr = VAR11 - 273.16
        VAR21_gr = VAR21 - 273.16
        VAR31_gr = VAR31 - 273.16
        VAR41_gr = VAR41 - 273.16
        
        VAR12_gr = VAR12 - 273.16
        VAR22_gr = VAR22 - 273.16
        VAR32_gr = VAR32 - 273.16
        VAR42_gr = VAR42 - 273.16
        
        lat = VAR11.lat
        lon = VAR11.lon
        
    elif var =='TAUX' or var == 'TAUY' or var == 'ICEFRAC' or var =='PRECT' or var =='NBP':
        VAR11_gr = VAR11
        VAR21_gr = VAR21
        VAR31_gr = VAR31
        VAR41_gr = VAR41
        
        VAR12_gr = VAR12
        VAR22_gr = VAR22
        VAR32_gr = VAR32
        VAR42_gr = VAR42
        
        lat = VAR11.lat
        lon = VAR11.lon
        
    else:
         VAR11_gr = regridder(VAR11)
         VAR11_gr=np.roll(VAR11_gr,-180)   
         VAR21_gr = regridder(VAR21)
         VAR21_gr=np.roll(VAR21_gr,-180) 
         VAR31_gr = regridder(VAR31)
         VAR31_gr=np.roll(VAR31_gr,-180)   
         VAR41_gr = regridder(VAR41)
         VAR41_gr=np.roll(VAR41_gr,-180) 
         
         VAR12_gr = regridder(VAR12)
         VAR12_gr=np.roll(VAR12_gr,-180)   
         VAR22_gr = regridder(VAR22)
         VAR22_gr=np.roll(VAR22_gr,-180) 
         VAR32_gr = regridder(VAR32)
         VAR32_gr=np.roll(VAR32_gr,-180)   
         VAR42_gr = regridder(VAR42)
         VAR42_gr=np.roll(VAR42_gr,-180) 
         
         lat = area_gr.lat
         lon = area_gr.lon
         
    #%% Plotting variables
    FS=20 # Base fontsize
        
    #%%
    # Directory for figures
    datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'
    
    # Select bounds for color scheme
    Vmn1 = ([-10,0,-0.2,0,0,-5e-9,-1e-9,0,32,0,0,7.6,280,280,0,0,0,0])
    Vmx1 = ([35,10,0.2,0.4,0.4,5e-9,1e-9,35,38,6,6e-7,8.15,360,380,0.25,1,1.2e-7,5e-8])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figures (some are extra)
    scen = '585'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR41_gr,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_585_HOS1.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR32_gr,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_585_CTL2.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR42_gr,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_585_HOS2.png', format='png', dpi=quality,bbox_inches='tight')
        
    #%%
    # Select bounds for color scheme
    Vmn1 = ([-10,0,-0.2,0,0,-5e-9,-1e-9,0,32,0,0,7.9,280,280,0,0,0,0])
    Vmx1 = ([35,10,0.2,0.4,0.4,5e-9,1e-9,35,38,6,6e-7,8.15,360,380,0.25,1,1.2e-7,5e-8])
    
    #%%
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figures (some are extra)
    scen = '126'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR21_gr,0,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_126_HOS1.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR12_gr,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_126_CTL2.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '126'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR22_gr,1,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_126_HOS2.png', format='png', dpi=quality,bbox_inches='tight')
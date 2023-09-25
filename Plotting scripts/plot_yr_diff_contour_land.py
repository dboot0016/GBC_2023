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

VARS= 'NBP','ALT','GPP','SR'
UNIT = '[kg C m$^{-2}$ s$^{-1}$]','[m]','[kg C m$^{-2}$ s$^{-1}$]','[kg C m$^{-2}$ s$^{-1}$]'
DESCR = 'NBP','Active layer thickness','Gross Primary Production','Soil Respiration'
CMAP = (['RdYlGn','copper_r',cm.cm.algae,'YlOrBr'])
FIG_NR = '9','S19','S18','S20'

lat1 = -90
lat2 = 90
lon1 = -179-60
lon2 = 179-60

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
def subplot(data1,data2,data3,data4,i,scen):
    if var == 'ALT':
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    
    else:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
        ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='50m',zorder=5)
    
    if var is not 'ALT':
        ax.add_feature(cfeature.OCEAN)

         
    if i == -1:
        im=plt.pcolormesh(lon,lat,data3-data1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference (CTL-' + str(scen)+')',fontsize=FS)
        
    elif i == -2:
        im=plt.pcolormesh(lon,lat,data4-data2,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference (HOS-' + str(scen)+')',fontsize=FS)

    else:    
        im=plt.pcolormesh(lon,lat,data4-data3,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference HOS-CTL (' + str(scen)+')',fontsize=FS)

    # Colorbar specifics
    cbar=plt.colorbar(im,orientation='horizontal',) 
    cbar.ax.set_xlabel('$\Delta$'+descr+' ' + unit, fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)

    
    if lon1 == (-179-60) and not var == 'ALT':
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

         
    #%% Plotting variables
    FS=20 # Base fontsize
        
    #%%
    # Directory for figures
    datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'
    
    # Select bounds for color scheme
    Vmn1=([-1.6e-9,-1,-3e-8,-9e-8])
    Vmx1=([1.6e-9,1,3e-8,9e-8])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
   
    # Plot figure
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR41_gr,VAR32_gr,VAR42_gr,8,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'f.png', format='png', dpi=quality,bbox_inches='tight')
    
    # Plot figure
    scen = '126'
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR21_gr,VAR12_gr,VAR22_gr,8,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'c.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    # Select bounds for color scheme
    Vmn1=([-4.6e-9,-1,-7e-8,-3e-8])
    Vmx1=([4.6e-9,1,7e-8,3e-8])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figure
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR41_gr,VAR32_gr,VAR42_gr,-1,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'e.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    # Select bounds for color scheme
    Vmn1=([-4.6e-9,-1,-2e-8,-9e-9])
    Vmx1=([4.6e-9,1,2e-8,9e-9])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figure
    scen = '126'
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR21_gr,VAR12_gr,VAR22_gr,-1,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'b.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    # Select bounds for color scheme
    Vmn1=([-4.6e-9,-1,-7e-8,-3e-8])
    Vmx1=([4.6e-9,1,7e-8,3e-8])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figure (this figure is not included in supplementary material)
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR41_gr,VAR32_gr,VAR42_gr,-2,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_585_diff.png', format='png', dpi=quality,bbox_inches='tight')
    
    #%%
    # Select bounds for color scheme
    Vmn1=([-4.6e-9,-1,-2e-8,-9e-9])
    Vmx1=([4.6e-9,1,2e-8,9e-9])
    
    vmn1=Vmn1[NN]
    vmx1=Vmx1[NN]
    
    # Plot figure  (this figure is not included in supplementary material)
    scen = '126'
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR21_gr,VAR12_gr,VAR22_gr,-2,scen)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'_extra_126_diff.png', format='png', dpi=quality,bbox_inches='tight')
    

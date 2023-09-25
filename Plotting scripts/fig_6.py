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

VARS= 'DIC',''
UNIT = '[mol m$^{-3}$]',''
DESCR = 'DIC',''
CMAP = (['RdYlBu_r',''])

NN = 0

scen = '585'   
reg = 'glob'

#%%
data1='/Users/daan/Documents/Snellius_stuff/Repository/Data/'   # Location of dataset(s) 
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
def MakeDataArray(data, z_t,lat, lon, dims):
        
        array = xr.DataArray(data = data,
                                 dims = ["z_t", "lat", "lon"],
                                 coords = dict(z_t = z_t,
                                               lat = lat,
                                               lon = lon));   
        return array
    
#%%
def subplot_diff(data1,data2,data3,data4,i,scen,basin):
    
    ax1 = ax[0]
    ax2 = ax[1]
    
    if i == -1:
        im=ax1.pcolormesh(lat,-d,data2-data1,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference '+str(basin)+' (CTL-' + str(scen)+')',fontsize=FS)
        ax1.contour(lat,-d,data2-data1,[0],colors='black')
        
        ax2.pcolormesh(lat,-d,data2-data1,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax2.contour(lat,-d,data2-data1,[0],colors='black')
        
    elif i == -2:
        im=ax1.pcolormesh(lat,-d,data4-data3,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference '+str(basin)+' (HOS-' + str(scen)+')',fontsize=FS)
        ax1.contour(lat,-d,data4-data3,[0],colors='black')
        
        ax2.pcolormesh(lat,-d,data4-data3,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax2.contour(lat,-d,data4-data3,[0],colors='black')
        
    else:    
        im=ax1.pcolormesh(lat,-d,data4-data2,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference HOS-CTL (' + str(scen)+'; '+str(basin)+')',fontsize=FS)
        ax1.contour(lat,-d,data4-data2,[0],colors='black')
        
        ax2.pcolormesh(lat,-d,data4-data2,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax2.contour(lat,-d,data4-data2,[0],colors='black')

    # Colorbar specifics
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar=fig.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel('$\Delta$'+descr+' ' + unit, fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    ax2.set_xlabel('Latitude [$^{\circ}$N]',fontsize=FS-2)
    fig.text(-0.04, 0.5, 'Depth [m]', va='center', rotation='vertical',fontsize=FS-2)
    #ax1.set_ylabel('Depth [m]',fontsize=FS-2)
    
    ax1.tick_params(axis='both', which='major', labelsize=FS-4)
    ax1.tick_params('x', labelbottom=False)
    ax2.tick_params(axis='both', which='major', labelsize=FS-4)

    ax1.set_ylim([-999,0])
    ax2.set_ylim([-max(d),-1000])
    
    ax1.set_xlim([-80,70])
    ax2.set_xlim([-80,70])
    
    ax1.grid()
    ax2.grid()

#%%    
def subplot_abs(data1,data2,i,exp,scen,basin):
    ax1 = ax[0]
    ax2 = ax[1]
    
    if i == 0:
        im=ax1.pcolormesh(lat,-d,data1,vmin=vmn1,vmax=vmx1,cmap=cmap)
        ax1.set_title('Average 2016-2020 ('+str(scen)+'-' + str(exp)+'; '+str(basin)+')',fontsize=FS)
        ax2.pcolormesh(lat,-d,data1,vmin=vmn1,vmax=vmx1,cmap=cmap)
        
    elif i == 1:
        im=ax1.pcolormesh(lat,-d,data2,vmin=vmn1,vmax=vmx1,cmap=cmap)
        ax1.set_title('Average 2096-2100 ('+str(scen)+'-' + str(exp)+'; '+str(basin)+')',fontsize=FS)
        ax2.pcolormesh(lat,-d,data2,vmin=vmn1,vmax=vmx1,cmap=cmap)
        
    # Colorbar specifics
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar=fig.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel('$\Delta$'+descr+' ' + unit, fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    ax2.set_xlabel('Latitude [$^{\circ}$N]',fontsize=FS-2)
    fig.text(-0.04, 0.5, 'Depth [m]', va='center', rotation='vertical',fontsize=FS-2)
    #ax1.set_ylabel('Depth [m]',fontsize=FS-2)
    
    ax1.tick_params(axis='both', which='major', labelsize=FS-4)
    ax1.tick_params('x', labelbottom=False)
    ax2.tick_params(axis='both', which='major', labelsize=FS-4)

    ax1.set_ylim([-999,0])
    ax2.set_ylim([-max(d),-1000])
    
    ax1.set_xlim([-80,70])
    ax2.set_xlim([-80,70])
    
    ax1.grid()
    ax2.grid()

#%%
var=VARS[NN]
unit=UNIT[NN]
descr=DESCR[NN]
cmap=CMAP[NN]

print(var)
      
#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_depth_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
VAR2_1=load_var2[var].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_depth_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0126.nc' )
VAR1_1=load_var1[var].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_depth_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
VAR4_1=load_var2[var].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_depth_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_0585.nc' )
VAR3_1=load_var1[var].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_depth_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
VAR2_2=load_var2[var].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_depth_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0126.nc' )
VAR1_2=load_var1[var].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_depth_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
VAR4_2=load_var2[var].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_depth_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_0585.nc' )
VAR3_2=load_var1[var].compute().squeeze()

d=load_var1['z_t']/100

#%%
VAR1_gr1 = regridder(VAR1_1)
VAR1_gr1=np.roll(VAR1_gr1,-180)   

VAR2_gr1 = regridder(VAR2_1)
VAR2_gr1=np.roll(VAR2_gr1,-180) 

VAR3_gr1 = regridder(VAR3_1)
VAR3_gr1=np.roll(VAR3_gr1,-180)   

VAR4_gr1 = regridder(VAR4_1)
VAR4_gr1=np.roll(VAR4_gr1,-180) 

VAR1_gr2 = regridder(VAR1_2)
VAR1_gr2=np.roll(VAR1_gr2,-180)   

VAR2_gr2 = regridder(VAR2_2)
VAR2_gr2=np.roll(VAR2_gr2,-180) 

VAR3_gr2 = regridder(VAR3_2)
VAR3_gr2=np.roll(VAR3_gr2,-180)   

VAR4_gr2 = regridder(VAR4_2)
VAR4_gr2=np.roll(VAR4_gr2,-180) 

lat = area_gr.lat
lon = area_gr.lon

#%%
read_var = MakeDataArray(data = VAR1_gr1, z_t=d, lat=lat, lon=lon, dims = "3d")
VAR1_gr1 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR2_gr1,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR2_gr1 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR3_gr1,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR3_gr1 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR4_gr1,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR4_gr1 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR1_gr2,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR1_gr2 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR2_gr2, z_t=d, lat=lat, lon=lon, dims = "3d")
VAR2_gr2 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR3_gr2,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR3_gr2 = read_var.to_dataset(name = 'DIC');

read_var = MakeDataArray(data = VAR4_gr2,  z_t=d, lat=lat, lon=lon, dims = "3d")
VAR4_gr2 = read_var.to_dataset(name = 'DIC');

#%%
VAR1_1, VAR2_1, VAR3_1, VAR4_1 = [],[],[],[]
VAR1_2, VAR2_2, VAR3_2, VAR4_2 = [],[],[],[]

#%%
RE = 6.371e6  # [m] Earth radius                           # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360    # grid size in x-direction

#%%
total_dx=((VAR1_gr1['DIC'].count(['lon']))*dx)

#%%
VAR1_glob1=(VAR1_gr1*dx).sum(['lon'])/total_dx
VAR2_glob1=(VAR2_gr1*dx).sum(['lon'])/total_dx
VAR3_glob1=(VAR3_gr1*dx).sum(['lon'])/total_dx
VAR4_glob1=(VAR4_gr1*dx).sum(['lon'])/total_dx

VAR1_glob2=(VAR1_gr2*dx).sum(['lon'])/total_dx
VAR2_glob2=(VAR2_gr2*dx).sum(['lon'])/total_dx
VAR3_glob2=(VAR3_gr2*dx).sum(['lon'])/total_dx
VAR4_glob2=(VAR4_gr2*dx).sum(['lon'])/total_dx

#%%
mask = (VAR1_gr1['DIC'][0,:,:])
mask1 = xr.zeros_like(mask)*mask
mask1[:,20:-110] = np.nan
mask1[70+90:,:] = np.nan
mask1[:95,-110:-70] = np.nan
mask1[110:140,:30] = np.nan
mask1[140:160,14:290] = np.nan
mask1[95:108,-110:-90] = np.nan
mask1[108:120,-110:-105] = np.nan
mask1[95:101,-110:-78] = np.nan

fig = plt.figure(figsize=(7, 5))
plt.contourf(mask1+1)

#%%
mask = (VAR1_gr1['DIC'][0,:,:])
mask2 = (xr.zeros_like(mask)*mask)
mask2 = mask2-(mask1.fillna(1))
mask2 = (mask2.where(mask2==-1))
mask2[66+90:,:]=np.nan
mask2[110:160,:34] = np.nan
mask2[130:160,:45] = np.nan
mask2[140:160,260:] = np.nan
mask2[96:110,278:] = np.nan

fig = plt.figure(figsize=(7, 5))
plt.contourf(mask2)
plt.colorbar()

#%%
mask1 = mask1+1
mask2 = mask2+2

#%%
total_dx=(((VAR1_gr1['DIC']*mask1).count(['lon']))*dx)

#%%
VAR1_atl1 = (VAR1_gr1*mask1*dx).sum(['lon'])/total_dx
VAR2_atl1 = (VAR2_gr1*mask1*dx).sum(['lon'])/total_dx
VAR3_atl1 = (VAR3_gr1*mask1*dx).sum(['lon'])/total_dx
VAR4_atl1 = (VAR4_gr1*mask1*dx).sum(['lon'])/total_dx

VAR1_atl2 = (VAR1_gr2*mask1*dx).sum(['lon'])/total_dx
VAR2_atl2 = (VAR2_gr2*mask1*dx).sum(['lon'])/total_dx
VAR3_atl2 = (VAR3_gr2*mask1*dx).sum(['lon'])/total_dx
VAR4_atl2 = (VAR4_gr2*mask1*dx).sum(['lon'])/total_dx

#%%
total_dx=(((VAR1_gr1['DIC']*mask2).count(['lon']))*dx)
VAR1_pac1 = (VAR1_gr1*mask2*dx).sum(['lon'])/total_dx
VAR2_pac1 = (VAR2_gr1*mask2*dx).sum(['lon'])/total_dx
VAR3_pac1 = (VAR3_gr1*mask2*dx).sum(['lon'])/total_dx
VAR4_pac1 = (VAR4_gr1*mask2*dx).sum(['lon'])/total_dx

VAR1_pac2 = (VAR1_gr1*mask2*dx).sum(['lon'])/total_dx
VAR2_pac2 = (VAR2_gr1*mask2*dx).sum(['lon'])/total_dx
VAR3_pac2 = (VAR3_gr1*mask2*dx).sum(['lon'])/total_dx
VAR4_pac2 = (VAR4_gr1*mask2*dx).sum(['lon'])/total_dx

#%% Plotting variables
FS=20
    
#%%
datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'


Basin = (['global','Atlantic','Indo-Pacific'])
for n in range(1):
    n = 1
    basin = Basin[n]
    
    if n == 0:
        VAR1_gr1 = VAR1_glob1['DIC']
        VAR2_gr1 = VAR2_glob1['DIC']
        VAR3_gr1 = VAR3_glob1['DIC']
        VAR4_gr1 = VAR4_glob1['DIC']
        
        VAR1_gr2 = VAR1_glob2['DIC']
        VAR2_gr2 = VAR2_glob2['DIC']
        VAR3_gr2 = VAR3_glob2['DIC']
        VAR4_gr2 = VAR4_glob2['DIC']
        
    elif n == 1:
        VAR1_gr1 = VAR1_atl1['DIC']
        VAR2_gr1 = VAR2_atl1['DIC']
        VAR3_gr1 = VAR3_atl1['DIC']
        VAR4_gr1 = VAR4_atl1['DIC']
        
        VAR1_gr2 = VAR1_atl2['DIC']
        VAR2_gr2 = VAR2_atl2['DIC']
        VAR3_gr2 = VAR3_atl2['DIC']
        VAR4_gr2 = VAR4_atl2['DIC']
        
    elif n == 2:
        VAR1_gr1 = VAR1_pac1['DIC']
        VAR2_gr1 = VAR2_pac1['DIC']
        VAR3_gr1 = VAR3_pac1['DIC']
        VAR4_gr1 = VAR4_pac1['DIC']
        
        VAR1_gr2 = VAR1_pac2['DIC']
        VAR2_gr2 = VAR2_pac2['DIC']
        VAR3_gr2 = VAR3_pac2['DIC']
        VAR4_gr2 = VAR4_pac2['DIC']
    
    
    Vmn1 = ([-0.16,-0.16,-0.16])
    
    vmn1=Vmn1[n]
    vmx1=-vmn1

    
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,8,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6f.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    scen = '126'
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr2,8,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6c.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    Vmn1 = ([-0.11,-0.11,-0.11])
    
    vmn1=Vmn1[n]
    vmx1=-vmn1
    
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,-1,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6e.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    scen = '126'
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr2,-1,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6b.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585'
    i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    
    subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,-2,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_diff_HOS_585.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr2,-2,scen,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_diff_HOS_126.png', format='png', dpi=quality,bbox_inches='tight')
    
    Vmn1 = ([2,2,2])
    Vmx1 = ([2.6,2.6,2.6])
    
    vmn1=Vmn1[n]
    vmx1=Vmx1[n]
    
    scen = '585'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR3_gr1,VAR3_gr2,0,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6d.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    scen = '585'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR4_gr1,VAR4_gr2,0,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_HOS_585_2016.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR3_gr1,VAR3_gr2,1,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_CTL_585_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR4_gr1,VAR4_gr2,1,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_HOS_585_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    
    Vmn1 = ([2,2,2])
    Vmx1 = ([2.6,2.6,2.6])
    
    vmn1=Vmn1[n]
    vmx1=Vmx1[n]
    
    scen = '126'
    exp = 'CTL'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR1_gr1,VAR1_gr2,0,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6a.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'HOS'
    i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR2_gr1,VAR2_gr2,0,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_HOS_126_2016.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'CTL'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR1_gr1,VAR1_gr2,1,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_CTL_126_2096.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126'
    exp = 'HOS'
    i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
    subplot_abs(VAR2_gr1,VAR2_gr2,1,scen,exp,basin)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'fig_6_extra_HOS_126_2096.png', format='png', dpi=quality,bbox_inches='tight')

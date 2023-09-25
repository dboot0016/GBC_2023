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

VARS= 'MOC',''
UNIT = '[Sv]',''
DESCR = 'AMOC',''
CMAP = ([cm.cm.curl,''])

NN = 0

scen = '585'   
reg = 'glob'
    
#%% Location data
data_snel1='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_126/'
data_snel2='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_126/'
data_snel3='/Users/daan/Documents/Snellius_stuff/Repository/Data/CTL_585/'
data_snel4='/Users/daan/Documents/Snellius_stuff/Repository/Data/HOS_585/' 

#%%
def subplot_diff(data1,data2,data3,data4,i,scen):
    
    ax1 = ax[0]
    ax2 = ax[1]
    
    if i == -1:
        im=ax1.pcolormesh(lat,-d,data2-data1,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference  (CTL-' + str(scen)+')',fontsize=FS)
        ax1.contour(lat,-d,data2-data1,[0],colors='black')
        
        ax2.pcolormesh(lat,-d,data2-data1,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax2.contour(lat,-d,data2-data1,[0],colors='black')
        
    elif i == -2:
        im=ax1.pcolormesh(lat,-d,data4-data3,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference  (HOS-' + str(scen)+')',fontsize=FS)
        ax1.contour(lat,-d,data4-data3,[0],colors='black')
        
        ax2.pcolormesh(lat,-d,data4-data3,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax2.contour(lat,-d,data4-data3,[0],colors='black')
        
    else:    
        im=ax1.pcolormesh(lat,-d,data4-data2,vmin=vmn1,vmax=vmx1,cmap='RdBu_r')
        ax1.set_title('Difference HOS-CTL (' + str(scen)+')',fontsize=FS)
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
    
    ax1.set_xlim([-34,70])
    ax2.set_xlim([-34,70])
    
    ax1.grid()
    ax2.grid()

#%%    
def subplot_abs(data1,data2,i,exp,scen):
    ax1 = ax[0]
    ax2 = ax[1]
    
    if i == 0:
        im=ax1.pcolormesh(lat,-d,data1,vmin=vmn1,vmax=vmx1,cmap=cmap)
        ax1.set_title('Average 2016-2020 ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
        ax2.pcolormesh(lat,-d,data1,vmin=vmn1,vmax=vmx1,cmap=cmap)
        
    elif i == 1:
        im=ax1.pcolormesh(lat,-d,data2,vmin=vmn1,vmax=vmx1,cmap=cmap)
        ax1.set_title('Average 2096-2100 ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
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
    
    ax1.set_xlim([-34,70])
    ax2.set_xlim([-34,70])
    
    ax1.grid()
    ax2.grid()
    

#%%
var=VARS[NN]
unit=UNIT[NN]
descr=DESCR[NN]
cmap=CMAP[NN]

print(var)
      
#%% Call on datasets (might need to change name of dataset)
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2_gr1=load_var2[var][:5,1,0,:,:].mean('time').compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1_gr1=load_var1[var][:5,1,0,:,:].mean('time').compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4_gr1=load_var2[var][:5,1,0,:,:].mean('time').compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3_gr1=load_var1[var][:5,1,0,:,:].mean('time').compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR2_gr2=load_var2[var][-5:,1,0,:,:].mean('time').compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_126)+'_0126.nc' )
VAR1_gr2=load_var1[var][-5:,1,0,:,:].mean('time').compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR4_gr2=load_var2[var][-5:,1,0,:,:].mean('time').compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end_585)+'_0585.nc' )
VAR3_gr2=load_var1[var][-5:,1,0,:,:].mean('time').compute().squeeze()

d=load_var1['moc_z']/100
lat=load_var1['lat_aux_grid']
     
#%% Plotting variables
FS=20
    
#%%
datadir = '/Users/daan/Documents/Snellius_stuff/Repository/Figures/'

vmn1=-10
vmx1=10

scen = '585'
i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,8,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22f.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr2,8,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22c.png', format='png', dpi=quality,bbox_inches='tight')

#%%
vmn1=-15
vmx1=15

scen = '585'
i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,-1,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22e.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr1,-1,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22b.png', format='png', dpi=quality,bbox_inches='tight')

scen = '585'
i = 8 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR3_gr1,VAR3_gr2,VAR4_gr1,VAR4_gr2,-2,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_diff_HOS_585.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_diff(VAR1_gr1,VAR1_gr2,VAR2_gr1,VAR2_gr2,-2,scen)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_diff_HOS_126.png', format='png', dpi=quality,bbox_inches='tight')

#%%
vmn1=-20
vmx1=20

scen = '585'
exp = 'CTL'
i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR3_gr1,VAR3_gr2,0,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22d.png', format='png', dpi=quality,bbox_inches='tight')

#%%
scen = '585'
exp = 'HOS'
i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR4_gr1,VAR4_gr2,0,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_HOS_585_2016.png', format='png', dpi=quality,bbox_inches='tight')

scen = '585'
exp = 'CTL'
i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR3_gr1,VAR3_gr2,1,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_CTL_585_2096.png', format='png', dpi=quality,bbox_inches='tight')

scen = '585'
exp = 'HOS'
i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR4_gr1,VAR4_gr2,1,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_HOS_585_2096.png', format='png', dpi=quality,bbox_inches='tight')

#%%
#vmn1 = ([-10,0,-0.2,-0.2,0,-5e-9,-1e-9,0,32,0,0,0,0,0,0,0,0,0,-30])
#vmx1 = ([35,10,0.2,0.2,0.4,5e-9,1e-9,35,38,6,6e-7,1e-7,0.25,0.06,0.12,0.005,0.06,0.005,30])

vmn1=-20
vmx1=20

scen = '126'
exp = 'CTL'
i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR1_gr1,VAR1_gr2,0,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22a.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
exp = 'HOS'
i = 0 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR2_gr1,VAR2_gr2,0,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_HOS_126_2016.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
exp = 'CTL'
i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR1_gr1,VAR1_gr2,1,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_CTL_126_2096.png', format='png', dpi=quality,bbox_inches='tight')

scen = '126'
exp = 'HOS'
i = 1 # Which decade; select -1 for Diff CTL last 5 yrs. w.r.t first 5 yrs -2 for DIFF HOS
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),gridspec_kw={'height_ratios': [1, 2],'hspace': 0.1})
subplot_abs(VAR2_gr1,VAR2_gr2,1,scen,exp)

if save_fig == 'yes':
    plt.savefig(datadir+'fig_S22_extra_HOS_126_2096.png', format='png', dpi=quality,bbox_inches='tight')


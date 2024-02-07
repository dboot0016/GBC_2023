import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats

#%%
def weighted_temporal_mean(ds, dt):
    month_length = ds.time.dt.days_in_month
    wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)
    
    obs = ds
    cond = obs.isnull()
    ones = xr.where(cond, 0.0, 1.0)

    obs_sum = (obs * wgts).resample(time=str(dt)+"AS").sum(dim="time")
    ones_out = (ones * wgts).resample(time=str(dt)+"AS").sum(dim="time")

    return obs_sum / ones_out

#%%
def trends(data1):
    a1 = np.zeros((16,))
    a1[0] = np.mean(data1[5:15],axis=0)-np.mean(data1[:5],axis=0)
    
    for i in range(15):
        a1[i+1] = np.mean(data1[i*10+15:i*10+25],axis=0)-np.mean(data1[i*10+5:i*10+15],axis=0)    
    
    return a1

#%%
def linear_regress(data1,data2,size):
    Data1 = data1.reshape(size,1)
    Data2 = data2
     
    model = LinearRegression().fit(Data1,Data2)
    r_sq = model.score(Data1,Data2)
    
    a = model.coef_
    b = model.intercept_
    
    params = np.append(model.intercept_,model.coef_)
    predictions = model.predict(Data1)
    new_X = np.append(np.ones((len(Data1),1)), Data1, axis=1)
    M_S_E = (sum((Data2-predictions)**2))/(len(new_X)-len(new_X[0]))
    v_b = M_S_E*(np.linalg.inv(np.dot(new_X.T,new_X)).diagonal())
    s_b = np.sqrt(v_b)
    t_b = params/ s_b
    p_val =[2*(1-stats.t.cdf(np.abs(i),(len(new_X)-len(new_X[0])))) for i in t_b]
    p_val = np.round(p_val,3)
    
    return p_val, r_sq, predictions

#%%
LW = 4
FS = 20
def plot(pval,rsq,prediction,data1,data2,descr1,descr2,unit1,unit2):
    y = np.arange(2026,2100,10)
    u1 = 'dec$^{-1}$'
        
    plt.plot(data1,prediction,linewidth=LW,color = 'black',label = 'R$^2$ = ' + str(np.round(rsq,2)) +'; p: ' +str(pval))
    
    im = plt.scatter(data1,data2,c=y,s=200,vmin = 2016,vmax=2100,cmap = 'inferno',label = '_nolegend_')
    cbar=plt.colorbar(im,orientation='vertical',) 
    cbar.ax.set_xlabel('Year', fontsize=16)
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    plt.xlabel(descr1 + ' ['+unit1+' '+u1+']',fontsize=FS-2)
    plt.ylabel(descr2 + ' ['+unit2+' '+u1+']',fontsize=FS-2)
    plt.title(str(model) + ': ' + str(mem),fontsize=FS)
    plt.xticks(fontsize=FS-4)
    plt.yticks(fontsize=FS-4)
    plt.grid()
    plt.legend(fontsize=FS-4)
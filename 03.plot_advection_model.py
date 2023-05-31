#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:43:43 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

# model = 'w0201'
path = '/Users/chingchen/Desktop/model/'
datapath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/figure/StagYY/'
mp4 = 1
labelsize = 30
bwith = 3
fig_advection_model = 1
fig_temperature_model = 0
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
header_list = ['r','Tmean','Tmin','Tmax','vrms','vmin','vmax',
               'vzabs','vzmin','vzmax','vhrms','vhmin','vhmax',
               'etalog','etamin','etamax','elog','emin','emax',
               'slog','smin','smax','whrms','whmin','whmax',
               'wzrms','wzmin','wzmax','drms','dmin','dmax',
               'enadv','endiff','enradh','enviscdiss','enadiabh',
               'cmean','cmin','cmax','rhomean','rhomin','rhomax',
               'airmean','airmin','airmax','primmean','primmin',
               'primmax','ccmean','ccmin','ccmax','fmeltmean',
               'fmeltmin','fmeltmax','metalmean','metalmin',
               'metalmax','gsmean','gsmin','gsmax','viscdisslog',
               'viscdissmin','viscdissmax','advtot','advdesc',
               'advasc','tcondmean','tcondmin','tcondmax']


model_list = ['w0201','w0204','w0207','w0210']
# model_list = ['w0202','w0205','w0208','w0211']
# model_list = ['w0203','w0206','w0209','w0212']
model_list = ['w0213','w0216','w0219','w0222']
# model_list =['w0210']
model_list = ['w0204','w0216','w0228']
# model_list = ['w0204','w0205','w0206']
model_list = ['w0207','w0208','w0209']
model_list = ['w0219','w0220','w0221']
model_list = ['w0801','w0802','w0803','w0804']
model_list = ['w0201','w0801','w0213','w0225']
# model_list = ['test03','test01','test04']
# label_list = ['Ea=1e5','Ea=1e6','Ea=1e7','Ea=1e8']
#label_list = ['f=0.7','f=0.5','f=0.3']
label_list = ['Ra = 1e5','Ra=3.2e4','Ra = 1e4','Ra = 1e3']
end_list = [1500,1560,1500,300]
if fig_advection_model:
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(12,12))
    for kk, model in enumerate(model_list):
        ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(end_list[kk])+'.txt',
                          sep = '\\s+',header = None,names = header_list)    
        ax.plot(ff.vzabs*ff.Tmean,ff.r,color = newcolors[kk],lw=5,label = label_list[kk])
        ax2.plot(ff.Tmean,ff.r,color = newcolors[kk],lw=5,label = label_list[kk])
        
    for aa in [ax,ax2]:
        aa.tick_params(labelsize=labelsize)
        aa.grid()
        aa.set_ylim(0,1)
        for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
    ax2.set_xlim(0,1)
    ax.set_xlabel('Temperature x Vz',fontsize = labelsize)
    ax2.set_xlabel('Temperature ',fontsize = labelsize)
    ax.set_ylabel('Depth',fontsize = labelsize)
    # ax2.legend(fontsize = 30)
    # ax.set_title('Snapshot = '+str(end),fontsize = labelsize)
        

end = 150
for kk, model in enumerate(model_list):
    lid_thickness = np.zeros(end)
    new_thickness = np.zeros(end)
    avg_temp = np.zeros(end)
    time = np.linspace(1,end,end)
    for i in range(1,end+1):
        ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(end_list[kk])+'.txt',
                          sep = '\\s+',header = None,names = header_list)    
        x = np.array(ff.vzabs*ff.Tmean)
        y = np.array(ff.r)
        smooth_d2 = np.gradient(np.gradient(x))
        infls = np.where(np.diff(np.sign(smooth_d2)))[0]
 
        if len(x[infls])==0:
            inf_point = 0
            lid_thickness[i-1]=0
            avg_temp[i-1] = 0
        else:
            inf_point = x[infls][-1]
            index_inf_point = [i for i, kk in enumerate(x) if kk == inf_point][0]
            a = (y[index_inf_point]-y[index_inf_point+1])/(x[index_inf_point]-x[index_inf_point+1])
            b = y[infls][-1]-a*inf_point
            line_x = np.linspace(0,300)
            line_y = a*line_x + b
            new_thickness[i-1] = line_y[0]
            avg_temp[i-1] = np.average(ff.Tmean[y<y[x==inf_point]])
            
    ax.plot(line_x,line_y,color =newcolors[kk],lw = 2 )
    ax.scatter(x[index_inf_point],y[index_inf_point],color = 'orange',s = 300)
    ax.axhline(y=line_y[0], color=newcolors[kk], linestyle='--',lw = 2)
    ax2.axhline(y=line_y[0], color=newcolors[kk], linestyle='--',lw = 2)  
    print(model,line_y[0])
    if fig_temperature_model:
        fig2,(ax3) = plt.subplots(1,1,figsize=(12,6))
        ax3.plot(time,1-new_thickness,color = newcolors[kk],lw=4,label = label_list[kk])
        for aa in [ax3]:
            aa.tick_params(labelsize=labelsize)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
            aa.grid()
        ax3.set_xlabel('Step',fontsize = labelsize)
        ax3.set_ylabel('Thickness',fontsize = labelsize)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 15:20:40 2023

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
fig_temperature_model = 1
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



model = 'w0201'

end = 1000
if fig_advection_model:
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(12,12))
    # for kk, model in enumerate(model_list):
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax.plot(ff.vzabs*ff.Tmean,ff.r,color = '#414F67',lw=5)
    ax2.plot(ff.Tmean,ff.r,color = '#414F67',lw=5)
        
    for aa in [ax,ax2]:
        aa.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
        aa.grid()
        aa.set_ylim(0,1)
    ax2.set_xlim(0,1)
    ax.set_xlabel('Temperature x Vz',fontsize = labelsize)
    ax2.set_xlabel('Temperature ',fontsize = labelsize)
    ax.set_ylabel('Depth',fontsize = labelsize)
    ax.set_title('Snapshot = '+str(end),fontsize = labelsize)
        

lid_thickness = np.zeros(end)
new_thickness = np.zeros(end)
avg_temp = np.zeros(end)
time = np.linspace(1,end,end)
for i in range(1,end+1):
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
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
        line_x = np.linspace(0,2400)
        line_y = a*line_x + b
        new_thickness[i-1] = line_y[0]
        avg_temp[i-1] = np.average(ff.Tmean[y<y[x==inf_point]])
        
ax.plot(line_x,line_y,color ='#4198B9',lw = 3 )
ax.scatter(x[index_inf_point],y[index_inf_point],color = 'orange',s = 200)
ax.set_xlim(0,500)
ax.axhline(y=line_y[0], color='r', linestyle='--',lw = 2)
ax2.axhline(y=line_y[0], color='r', linestyle='--',lw = 2)

if fig_temperature_model:
    fig2,(ax3) = plt.subplots(1,1,figsize=(12,6))
    ax3.plot(time,1-new_thickness,color = '#4198B9',lw=4)
    for aa in [ax3]:
        aa.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
        aa.grid()
    ax3.set_xlabel('Step',fontsize = labelsize)
    ax3.set_ylabel('Thickness',fontsize = labelsize)
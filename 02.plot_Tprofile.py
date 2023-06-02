#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:01:43 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    modelpath = '/Users/chingchen/Desktop/model/'
else:
    path = '/lfs/jiching/data/'
    workpath = '/lfs/jiching/ScalingLaw_model/'
    modelpath = '/lfs/jiching/ScalingLaw_model/'
    figpath = '/lfs/jiching/figure/'

labelsize = 30
bwith = 3

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

model_information = pd.read_csv(workpath+'average_data.csv',sep=',') 
for jj in range(len(model_information)):
    model = model_information.model[jj]
    time_window2 = model_information.time_window2[jj]
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(time_window2*100)+'.txt',
                          sep = '\\s+',header = None,names = header_list)   
    x = np.array(ff.vzabs*ff.Tmean)
    y = np.array(ff.r)
    smooth_d2 = np.gradient(np.gradient(x))
    infls = np.where(np.diff(np.sign(smooth_d2)))[0]
    if len(x[infls])==0:
            inf_point = 0
            lid_thickness=0
            avg_temp = 0
    else:
        inf_point = x[infls][-1]
        index_inf_point = [i for i, kk in enumerate(x) if kk == inf_point][0]
        a = (y[index_inf_point]-y[index_inf_point+1])/(x[index_inf_point]-x[index_inf_point+1])
        b = y[infls][-1]-a*inf_point
        line_x = np.linspace(0,300)
        line_y = a*line_x + b
        new_thickness = line_y[0]
        avg_temp= np.average(ff.Tmean[y<y[x==inf_point]])
      
#    ax.plot(line_x,line_y,color =newcolors[kk],lw = 2 )
#    ax.scatter(x[index_inf_point],y[index_inf_point],color = 'orange',s = 300)
#    ax.axhline(y=line_y[0], color=newcolors[kk], linestyle='--',lw = 2)
#    ax2.axhline(y=line_y[0], color=newcolors[kk], linestyle='--',lw = 2)  
    print(model,line_y[0])

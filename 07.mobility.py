#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 23:13:08 2023

@author: chingchen
"""

import math
import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"


figpath = '/Users/chingchen/Desktop/figure/StagYY/'
datapath = '/Users/chingchen/Desktop/data/'
path = '/Users/chingchen/Desktop/model/'
model = 'w1009'
model_list = ['w0204','w0802','w0216']#,'w0228']
model_list = ['w1009']
# model_list = ['w0219','w0220','w0221']
# model_list = ['w0201','w0204','w0207']
# model_list = ['w0213','w0216','w0219','w0222']
# model_list = ['w0801','w0802','w0803','w0804']
# model_list = ['w0201','w0801','w0213']#,'w0225']
label_list = ['Ra = 1e5','Ra = 3.2e4','Ra = 1e4','Ra = 1e3']
labelsize = 30
bwith = 3
fig_mobility = 1
surface_velocity=0

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
header_list_timedat = ['istep','time','F_top','F_bot','Tmin',
               'Tmean','Tmax','Vmin','Vrms','Vmax','eta_min',
               'eta_mean','eta_max','ra_eff','Nu_top','Nu_bot',
               'C_min','C_mean','C_max','F_mean','F_max',
               'erupt_rate','erupta','erupt_heatflux',
               'entrainment','Cmass_error','H_int',
               'r_innercore','Tsurf','Tcmb']
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

end = 1500
if fig_mobility:
    fig4,(ax4) = plt.subplots(1,1,figsize=(12,8))
    kk = 3
    for kk, model in enumerate(model_list):
        file = datapath+model+'_time.dat'
        ff = pd.read_csv(file,sep = '\\s+')  
        time_vrms = np.average(ff.Vrms.tolist()[-500:])
        vsurf = np.zeros(end-1)
        vrms = np.zeros(end-1)
        snaptime = np.zeros(end-1)
        for uu in range(1,end):
            ffsnap = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(uu)+'.txt',
                          sep = '\\s+',header = None,names = header_list)
            vsurf[uu-1] = ffsnap.vhrms[127]
            vrms[uu-1] = np.average(ffsnap.vrms)
            snaptime[uu-1] = uu/1000
        ax4.plot(snaptime,vsurf/time_vrms,color = newcolors[kk],lw=5,label=label_list[kk]) #time.dat
        # ax4.plot(snaptime,vsurf/vrms,color = newcolors[kk],lw=5,label=model) #rprof
    ax4.set_ylim(0.0,0.01)
    ax4.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
        ax4.spines[axis].set_linewidth(bwith)
    ax4.grid()
    ax4.set_xlabel('time',fontsize = labelsize)
    ax4.legend(fontsize = 20)
    ax4.set_ylabel('Mobility',fontsize=labelsize)
if surface_velocity:
    fig5,(ax6) = plt.subplots(1,1,figsize=(12,4))
    vsurf = np.zeros(end-1)
    for uu in range(1,end):
        ffsnap = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(uu)+'.txt',
                      sep = '\\s+',header = None,names = header_list)
        vsurf[uu-1] = ffsnap.vhrms[127]
        snaptime[uu-1] = uu/1000
    ax6.plot(snaptime,vsurf,color = newcolors[kk],lw=5,label=model)
    # ax6.set_ylim(0.9,1)
    # ax4.set_xlim(0,10)
    ax6.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
        ax6.spines[axis].set_linewidth(bwith)
    ax6.grid()

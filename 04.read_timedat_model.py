#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:02:26 2023

@author: chingchen
"""
import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

# model = 'w0201'
path = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/StagYY/'
mp4 = 1
labelsize = 30
bwith = 3
fig_Nu_t = 0
fig_F_t  = 0
fig_T    = 1

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
header_list = ['istep','time','F_top','F_bot','Tmin',
               'Tmean','Tmax','Vmin','Vrms','Vmax','eta_min',
               'eta_mean','eta_max','ra_eff','Nu_top','Nu_bot',
               'C_min','C_mean','C_max','F_mean','F_max',
               'erupt_rate','erupta','erupt_heatflux',
               'entrainment','Cmass_error','H_int',
               'r_innercore','Tsurf','Tcmb']

#model_list = ['w0201','w0204','w0207','w0210']
# model_list = ['w0202','w0205','w0208','w0211']
model_list = ['w0207','w0208','w0209']
# model_list = ['w0201','w0213','w0225']
#model_list = ['w0213','w0216','w0219','w0222']
#model_list = ['w0204','w0216','w0228']
# model_list = ['w0214','w0217','w0220','w0223']
# model_list = ['w0203','w0206','w0209','w0212']
model_list = ['w0219','w0220','w0221']
label_list = ['Ea=1e5','Ea=1e6','Ea=1e7','Ea=1e8']
label_list = ['f=0.7','f=0.5','f=0.3']
#label_list = ['Ra=10$^5$','Ra=10$^4$','Ra=10$^3$']
path = '/Users/chingchen/Desktop/data/'

# for kk, model in enumerate(model_list):
if fig_Nu_t:
    fig,(ax) = plt.subplots(1,1,figsize=(12,6))
    for kk, model in enumerate(model_list):
        file = path+model+'_time.dat'
        ff = pd.read_csv(file,sep = '\\s+')    
        ax.plot(ff.time, ff.Nu_top,color = newcolors[kk],lw=5,label = 'Nu_top')
        ax.plot(ff.time, ff.Nu_bot,color = newcolors[kk+1],lw=3,label = 'Nu_bot')
    
    ax.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(bwith)
    ax.grid()
    # ax.set_xlim(0.09,0.12)
    # ax.set_ylim(7,8.2)
    ax.set_xlabel('t',fontsize = labelsize)
    ax.set_ylabel('heat flux',fontsize = labelsize)
    # ax.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
if fig_F_t:
    fig2,(ax2) = plt.subplots(1,1,figsize=(12,6))
    f = 0.5
    for kk, model in enumerate(model_list):
        file = path+model+'_time.dat'
        ff = pd.read_csv(file,sep = '\\s+')    
        # ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
                          # sep = '\\s+',header = None,names = header_list)    
        ax2.plot(ff.time,ff.F_bot,color = newcolors[kk],lw=3,label = label_list[kk])
        # ax2.plot(ff.Tmean,ff.r,color = newcolors[kk],lw=5,label = label_list[kk])
    
       
    # ax2.plot(ff.time, ff.F_bot*(f**2),color = newcolors[kk+1],lw=3,label = 'F_bot')
    # ax2.plot(ff.time, ff.F_top,color = newcolors[kk],lw=5,label = 'F_top')
        
    ax2.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
                ax2.spines[axis].set_linewidth(bwith)
    ax2.grid()
    # ax.set_xlim(0.09,0.12)
    # ax2.set_ylim(0,6)
    ax2.set_xlabel('t',fontsize = labelsize)
    ax2.set_ylabel('heat flux',fontsize = labelsize)
    ax2.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
if fig_T :
    fig3,(ax3) = plt.subplots(1,1,figsize=(12,6))
    for kk, model in enumerate(model_list):
        file = path+model+'_time.dat'
        ff = pd.read_csv(file,sep = '\\s+')
        ax3.plot(ff.time,ff.Tmean,color = newcolors[kk],lw=3,label = label_list[kk])
    # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk],lw=5,label = 'Tmean')
    # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk+1],lw=3,label = 'F_bot')
    
    ax3.tick_params( labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
                ax3.spines[axis].set_linewidth(bwith)
    ax3.grid()
    ax3.set_xlim(0,0.6)
    # ax.set_ylim(7,8.2)
    ax3.set_xlabel('time',fontsize = labelsize)
    ax3.set_ylabel('Temperature',fontsize = labelsize)
    ax3.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
# for qq in range(1000,1050):
    # print(ff.time[qq]-ff.time[qq-1])
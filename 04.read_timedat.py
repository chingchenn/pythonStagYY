#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 09:43:59 2023

@author: chingchen
"""
import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"


figpath = '/Users/chingchen/Desktop/figure/StagYY/'
mp4 = 1
labelsize = 30
bwith = 3
fig_Nu_t = 0
fig_F_t  = 0
fig_T    = 1
fig_mobility = 0

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

path = '/Users/chingchen/Desktop/data/'
model = 'h01'
f = 0.7
if model == 'w0213':
    xmin_Nu, xmax_Nu = 0.5,1.5
    ymin_Nu, ymax_Nu = 2,4
    ymin_F,ymax_F = 1,4
    f = 0.7
if model == 'w0216':
    xmin_Nu, xmax_Nu = 0,1.5
    ymin_Nu, ymax_Nu = 2,3
    ymin_F,ymax_F = 1,3
    f = 0.7
if model == 'w0201':
    xmin_Nu, xmax_Nu = 1.3,1.4
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,4
    f = 0.7
if model == 'w0204':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,4
    f = 0.7
if model == 'w0207':
    xmin_Nu, xmax_Nu = 0.0,0.3
    ymin_Nu, ymax_Nu = 4,10
    ymin_F,ymax_F = 3,8
    f = 0.7
if model == 'w0210':
    xmin_Nu, xmax_Nu = 0.05,0.2
    ymin_Nu, ymax_Nu = 4,10
    ymin_F,ymax_F = 3,10
    f = 0.7
elif model =='w0209':
    xmin_Nu, xmax_Nu = 0,0.3
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 1,5
    f = 0.3
    xmin_Nu, xmax_Nu = 0.26,0.29
    ymin_Nu, ymax_Nu = 4.9,5.1
    ymin_F,ymax_F = 1.49,1.52
if model == 'w0219':
    xmin_Nu, xmax_Nu = 0.0,1
    ymin_Nu, ymax_Nu = 2,6
    ymin_F,ymax_F = 0,4
    f = 0.7
if model == 'w0220':
    xmin_Nu, xmax_Nu = 0.0,1
    ymin_Nu, ymax_Nu = 2,6
    ymin_F,ymax_F = 0,4
    f = 0.5
if model == 'w0221':
    xmin_Nu, xmax_Nu = 0.0,1
    ymin_Nu, ymax_Nu = 2,6
    ymin_F,ymax_F = 0,4
    f = 0.3
elif model =='w0228':
    xmin_Nu, xmax_Nu = 0,0.3
    ymin_Nu, ymax_Nu = 0,4
    ymin_F,ymax_F = 0,2
    f = 0.7

xmin_Nu, xmax_Nu = 0,2.5
ymin_Nu, ymax_Nu = 1,10
ymin_F,ymax_F = 2,10

file = path+model+'_time.dat'
ff = pd.read_csv(file,sep = '\\s+')    
# for kk, model in enumerate(model_list):
if fig_Nu_t:
    fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
    kk = 1
    ax.plot(ff.time, ff.Nu_bot,color = newcolors[kk+1],lw=3,label = 'Nu_bot')   
    ax.plot(ff.time, ff.Nu_top,color = newcolors[kk],lw=3,label = 'Nu_top')
        
    ax.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(bwith)
    ax.grid()
    ax.set_xlim(xmin_Nu,xmax_Nu)
    ax.set_ylim(ymin_Nu, ymax_Nu)
    # ax.set_xlabel('time',fontsize = labelsize)
    ax.set_ylabel('heat flux',fontsize = labelsize)
    ax.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
# if fig_F_t:
    # fig2,(ax2) = plt.subplots(1,1,figsize=(12,6))
    kk = 1
    
       
    ax2.plot(ff.time, ff.F_bot*(f**2),color = newcolors[kk+1],lw=3,label = 'F_bot')
    ax2.plot(ff.time, ff.F_top,color = newcolors[kk],lw=3,label = 'F_top')
        
    ax2.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(bwith)
    ax2.grid()
    ax2.set_xlim(xmin_Nu,xmax_Nu)
    ax2.set_ylim(ymin_F,ymax_F)
    ax2.set_xlabel('time',fontsize = labelsize)
    ax2.set_ylabel('heat flux',fontsize = labelsize)
    ax2.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
if fig_T :
    fig3,(ax3) = plt.subplots(1,1,figsize=(12,6))
    kk = 1
       
    ax3.plot(ff.time, ff.Tmean,color = newcolors[kk],lw=5,label = 'Tmean')
    # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk+1],lw=3,label = 'F_bot')
    
    ax3.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
        ax3.spines[axis].set_linewidth(bwith)
    ax3.grid()
    # ax.set_xlim(0.09,0.12)
    ax3.set_ylim(0.78,0.828)
    ax3.set_xlabel('time',fontsize = labelsize)
    ax3.set_ylabel('Temperature',fontsize = labelsize)
    ax3.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
for qq in range(1000,1002):
    print(ff.time[qq]-ff.time[qq-1])
if fig_mobility:
    fig4,(ax4) = plt.subplots(1,1,figsize=(12,6))
    kk = 1
       
    ax4.plot(ff.time, ff.Tmean,color = newcolors[kk],lw=5,label = 'Tmean')
    # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk+1],lw=3,label = 'F_bot')
    
    ax4.tick_params(labelsize=labelsize)
    for axis in ['top','bottom','left','right']:
        ax4.spines[axis].set_linewidth(bwith)
    ax4.grid()
    # ax.set_xlim(0.09,0.12)
    # ax.set_ylim(7,8.2)
    ax4.set_xlabel('time',fontsize = labelsize)
    ax4.set_ylabel('Temperature',fontsize = labelsize)
    ax4.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 15:31:53 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"


figpath = '/Users/chingchen/Desktop/figure/StagYY/'
labelsize = 20
bwith = 3
fig_Nu_t = 1
fig_T    = 0
fig_error = 0

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
model = 'w1012'


model_information = pd.read_csv('/Users/chingchen/Desktop/StagYY_Works/model_information.csv',sep=',') 
f = np.array(model_information.f[model_information.model==model])[0]
Ea = np.array(model_information.Ea[model_information.model==model])[0]
Ra0 = np.array(model_information.Ra0[model_information.model==model])[0]
    
xmin_Nu, xmax_Nu = 0,2.5
ymin_Nu, ymax_Nu = 1,7
ymin_F,ymax_F = 1,7

# ymin_Nu, ymax_Nu = 2.5,3.2
# ymin_F,ymax_F = 2,2.6
time_window1,time_window2 = 2,2.1
error_y1, error_y2 = -0.3,0.3

if model == 'w0201':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,4
    f = 0.7
    time_window1,time_window2 = 2,2.1
    error_y1, error_y2 = -0.3,0.3
if model == 'w0204':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,4
    f = 0.7
    time_window1,time_window2 = 1.2,1.3
    error_y1, error_y2 = -0.6,0.6
if model == 'w0210':
    xmin_Nu, xmax_Nu = 0.05,2.5
    ymin_Nu, ymax_Nu = 4,10
    ymin_F,ymax_F = 3,10
    f = 0.7
    time_window1,time_window2 = 2,2.1
    error_y1, error_y2 = -0.3,0.3
if model == 'w0213':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 2,4
    ymin_F,ymax_F = 1,4
    f = 0.7
    time_window1,time_window2 = 2,2.1
    error_y1, error_y2 = -0.3,0.3
if model == 'w0219':
    xmin_Nu, xmax_Nu = 0.0,2.5
    ymin_Nu, ymax_Nu = 2,6
    ymin_F,ymax_F = 0,4
    f = 0.7
    time_window1,time_window2 = 1.5,1.7
    error_y1, error_y2 = -0.2,0.2
if model == 'w0222':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 0,10
    ymin_F,ymax_F = 0,4
    f = 0.7
    time_window1,time_window2 = 1,1.2
    error_y1, error_y2 = -0.3,0.3
if model == 'w0225':
    xmin_Nu, xmax_Nu = 0,2.5
    ymin_Nu, ymax_Nu = 0,10
    ymin_F,ymax_F = 0,4
    f = 0.7
    time_window1,time_window2 = 2,2.1
    error_y1, error_y2 = -0.3,0.3
if model == 'w0801':
    xmin_Nu, xmax_Nu =0,2.5
    ymin_Nu, ymax_Nu = 2,4
    ymin_F,ymax_F = 1,4
    f = 0.7
    time_window1,time_window2 = 2,2.1
    error_y1, error_y2 = -0.1,0.1
if model == 'w0802':
    xmin_Nu, xmax_Nu =0,2.5
    ymin_Nu, ymax_Nu = 2,4
    ymin_F,ymax_F = 1,4
    f = 0.7
    time_window1,time_window2 = 1.5,1.7
    error_y1, error_y2 = -0.3,0.3
if model == 'w0804':
    xmin_Nu, xmax_Nu =0,2.5
    ymin_Nu, ymax_Nu = 2,4
    ymin_F,ymax_F = 1,4
    f = 0.7
    time_window1,time_window2 = 1.5,1.7
    error_y1, error_y2 = -0.3,0.3
if model == 'w0809':
    xmin_Nu, xmax_Nu =0,2.5
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,5
    f = 0.8
    time_window1,time_window2 = 1.5,1.7
    error_y1, error_y2 = -0.3,0.3
if model == 'w0810':
    xmin_Nu, xmax_Nu =0,2.5
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 3,5
    f = 0.8
    time_window1,time_window2 = 0.4,0.5
    error_y1, error_y2 = -0.3,0.3    
if model == 'w0216':
    xmin_Nu, xmax_Nu = 0,1.5
    ymin_Nu, ymax_Nu = 2,3
    ymin_F,ymax_F = 1,3
    f = 0.7
if model == 'w0207':
    xmin_Nu, xmax_Nu = 0.0,0.3
    ymin_Nu, ymax_Nu = 4,10
    ymin_F,ymax_F = 3,8
    f = 0.7
if model =='w0209':
    xmin_Nu, xmax_Nu = 0,0.3
    ymin_Nu, ymax_Nu = 4,6
    ymin_F,ymax_F = 1,5
    f = 0.3
    xmin_Nu, xmax_Nu = 0.26,0.29
    ymin_Nu, ymax_Nu = 4.9,5.1
    ymin_F,ymax_F = 1.49,1.52
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
if model =='w0228':
    xmin_Nu, xmax_Nu = 0,0.3
    ymin_Nu, ymax_Nu = 0,4
    ymin_F,ymax_F = 0,2
    f = 0.7
file = path+model+'_time.dat'
ff = pd.read_csv(file,sep = '\\s+')    
tt1=[]
tt2=[]
anb=[]
ant=[]
afb=[]
aft=[]
for tt in np.arange(0.0,2.5,0.01):
    time_window1=tt
    time_window2=tt+0.1
    if len(ff.F_bot[(ff.time>time_window1) * (ff.time<time_window2)])==0:
        break
    average_Nu_bot = np.median(ff.Nu_bot[(ff.time>time_window1) * (ff.time<time_window2)])
    average_Nu_top = np.median(ff.Nu_top[(ff.time>time_window1) * (ff.time<time_window2)])
    average_fl_bot = np.median(ff.F_bot[(ff.time>time_window1) * (ff.time<time_window2)])
    average_fl_top = np.median(ff.F_top[(ff.time>time_window1) * (ff.time<time_window2)])
    
    tt1.append(np.round(time_window1,2))
    tt2.append(np.round(time_window1,2))
    anb.append(average_Nu_bot)
    ant.append(average_Nu_top)
    afb.append(average_fl_bot)
    aft.append(average_fl_top)

tt1 = np.array(tt1)
tt2 = np.array(tt2)
fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
ax.scatter(tt1,anb,color=newcolors[5])
ax2.scatter(tt1,afb,color=newcolors[5])

time_window1=np.max(tt1[(anb>np.mean(anb)+np.std(anb))])+0.05
time_window2=np.max(tt2[(anb<np.mean(anb)+np.std(anb))])-0.05
for aa in [ax,ax2]:
    aa.tick_params(labelsize=labelsize)
    aa.set_xlim(0,2.5)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
ax2.set_ylim(5,5.12)
ax.set_ylim(4,4.12)
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
if fig_error:
    kk = 1
    fig2,(ax5,ax6) = plt.subplots(2,1,figsize=(12,10))
    error_Nu_bot = ff.Nu_bot-average_Nu_bot
    error_Nu_top = ff.Nu_top-average_Nu_top
    error_fl_bot = ff.F_bot-average_fl_bot
    error_fl_top = ff.F_top-average_fl_top
    ax5.plot(ff.time, error_Nu_bot,color = newcolors[kk+1],lw=3,label = 'Nu_bot')   
    ax5.plot(ff.time, error_Nu_top,color = newcolors[kk],lw=3,label = 'Nu_top')
        
    ax5.tick_params(axis='x', labelsize=labelsize)
    ax5.tick_params(axis='y', labelsize=labelsize)
    ax5.spines['bottom'].set_linewidth(bwith)
    ax5.spines['top'].set_linewidth(bwith)
    ax5.spines['right'].set_linewidth(bwith)
    ax5.spines['left'].set_linewidth(bwith)
    ax5.grid()
    # ax5.set_xlim(xmin_Nu,xmax_Nu)
    # ax5.set_ylim(ymin_Nu, ymax_Nu)
    ax5.set_xlim(xmin_Nu,xmax_Nu)
    ax5.set_ylim(error_y1, error_y2)
    # ax.set_xlabel('time',fontsize = labelsize)
    ax5.set_ylabel('heat flux',fontsize = labelsize)
    ax5.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
     
    # ax6.plot(ff.time, ff.F_bot*(f**2),color = newcolors[kk+1],lw=3,label = 'F_bot')
    # ax6.plot(ff.time, ff.F_top,color = newcolors[kk],lw=3,label = 'F_top')
    
    ax6.plot(ff.time, error_fl_bot*(f**2),color = newcolors[kk+1],lw=3,label = 'F_bot')
    ax6.plot(ff.time, error_fl_top,color = newcolors[kk],lw=3,label = 'F_top')
    
    ax6.tick_params(axis='x', labelsize=labelsize)
    ax6.tick_params(axis='y', labelsize=labelsize)
    ax6.spines['bottom'].set_linewidth(bwith)
    ax6.spines['top'].set_linewidth(bwith)
    ax6.spines['right'].set_linewidth(bwith)
    ax6.spines['left'].set_linewidth(bwith)
    ax6.grid()
    ax6.set_ylim(error_y1, error_y2)
    ax6.set_xlim(xmin_Nu,xmax_Nu)
    # ax6.set_ylim(ymin_F,ymax_F)
    ax6.set_xlabel('time',fontsize = labelsize)
    ax6.set_ylabel('heat flux',fontsize = labelsize)

    
if fig_T :
    fig3,(ax3) = plt.subplots(1,1,figsize=(12,6))
    kk = 1
       
    ax3.plot(ff.time, ff.Tmean,color = newcolors[kk],lw=5,label = 'Tmean')
    # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk+1],lw=3,label = 'F_bot')
    
    ax3.tick_params(axis='x', labelsize=labelsize)
    ax3.tick_params(axis='y', labelsize=labelsize)
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax3.grid()
    ax3.set_xlim(0,2.5)
    ax3.set_ylim(0.5,1)
    ax3.set_xlabel('time',fontsize = labelsize)
    ax3.set_ylabel('Temperature',fontsize = labelsize)
    ax3.legend(fontsize = 20)
    # ax.set_title('t = '+str(i),fontsize = labelsize)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:41:20 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import time
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"


figpath = '/Users/chingchen/Desktop/figure/StagYY/'
labelsize = 20
bwith = 3


path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/'
modelpath = '/Users/chingchen/Desktop/model/'
model = 'w0816_192'

normal = 1
internal = 0

if normal:
    model_information = pd.read_csv(workpath+'model_information.csv',sep=',')
if internal:
    model_information = pd.read_csv(workpath+'model_internal_heating.csv',sep=',')  
#model_information = pd.read_csv('/Users/chingchen/Desktop/StagYY_Works/model_internal_heating.csv',sep=',')
fig_ftime = 1
choosen_time_window = 0
plot_average = 1
save_data = 0
rprof_header_list = ['r','Tmean','Tmin','Tmax','vrms','vmin','vmax',
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
for jj in range(len(model_information)):
    if model_information.model[jj]== model:
        break 
start = time.time()
f = model_information.f[jj]
Ea = model_information.Ea[jj]
Ra0 = model_information.Ra0[jj]


if internal:
    H = model_information.H[jj]

file = path+model+'_time.dat'
ff = pd.read_csv(file,sep = '\\s+')    
kk = np.arange(0.1,2.5,0.01) # time shift = 0.02
tt1=np.zeros(len(kk))
tt2=np.zeros(len(kk))
anb=np.zeros(len(kk))
ant=np.zeros(len(kk))
afb=np.zeros(len(kk))
aft=np.zeros(len(kk))

H = 1.0
#f=0
if fig_ftime:
    fig,(ax2) = plt.subplots(1,1,figsize=(12,8))
    if normal:
        ax2.plot(ff.time,ff.F_bot*(f**2),color='#AE6378',label = 'F_bot')
    if internal:
        Ftop_pred = f**2 * ff.F_bot + (1+f+f**2)/3 *H
        Ftop_pred = ff.F_bot + H
        ax2.plot(ff.time,Ftop_pred,color='#AE6378',label = 'F_bot')
    ax2.plot(ff.time,ff.F_top,color='#35838D',label = 'F_top',lw=3)
    #ax2.scatter(tt1,afb,color='#35838D')
    ax2.set_title(model,fontsize=labelsize)
    ax2.set_ylabel('F$_{bot}$',fontsize = labelsize)
    for aa in [ax2]:
        aa.legend(fontsize =labelsize)
        aa.tick_params(labelsize=labelsize)
        ax2.set_ylim(2,6)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)

# calculate the average Nu and heat flow in given time window    
time_window1=0.02
time_window2 = 0.07
average_Nu_bot = np.average(ff.Nu_bot[(ff.time>time_window1)&(ff.time<time_window2)])
average_Nu_top = np.average(ff.Nu_top[(ff.time>time_window1)&(ff.time<time_window2)])
average_fl_bot = np.average(ff.F_bot[(ff.time>time_window1)&(ff.time<time_window2)])
average_fl_top = np.average(ff.F_top[(ff.time>time_window1)&(ff.time<time_window2)])


print("total time taken this loop: ", time.time() - start)
### Find Average Temperature and the Stagnant Lid Thickness

end = int((time_window2-0.001)*1000)



lid_thickness = np.zeros(end)
new_thickness = np.zeros(end)
avg_temp = np.zeros(end)
time = np.linspace(1,end,end)

for i in range(1,end+1):
    ffrof = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = rprof_header_list)    
    x = np.array(ffrof.vzabs*ffrof.Tmean)
    y = np.array(ffrof.r)
    smooth_d2 = np.gradient(np.gradient(x))
    infls = np.where(np.diff(np.sign(smooth_d2)))[0]
 
    if len(x[infls])!=0:
        inf_point = x[infls][-1]
        index_inf_point = [i for i, kk in enumerate(x) if kk == inf_point][0]
        a = (y[index_inf_point]-y[index_inf_point+1])/(x[index_inf_point]-x[index_inf_point+1])
        b = y[infls][-1]-a*inf_point
        line_x = np.linspace(0,4000)
        line_y = a*line_x + b
        new_thickness[i-1] = line_y[0]
        avg_temp[i-1] = np.average(ffrof.Tmean[y<y[x==inf_point]])



rprof_ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
              sep = '\\s+',header = None,names = rprof_header_list) 
Tmm = np.average(rprof_ff.Tmean[(rprof_ff.r>0.3)*(rprof_ff.r<0.4)])
gamma = np.log(np.array(Ea))
rsurf = Ra0/np.exp(gamma/2)
raeff = rsurf*np.exp(gamma*Tmm) 


if internal:
    dlid = 1-np.average(new_thickness[50:])
    Ftop = average_fl_top
    Fbot = average_fl_bot
    
    
    flid = 1-(1-f)*dlid
    Tlid = dlid/flid*(Ftop-1/6*H*(2-flid-flid**2)/(1-f))
    Ur = (1+f+f**2)/3*H/Ftop
print('##############################################################')
if normal:
    print(model,round(rsurf,3),round(gamma,3),f,round(Tmm,3),round(average_fl_top,3),round(average_fl_bot,3),format(raeff,'.3e'),round(time_window1,3),round(time_window2,3),format(Ea,'.1e'),)
if internal:
    print(model,round(rsurf,3),f,format(Ea,'.1e'),H,round(Tmm,3),round(average_fl_top,3),round(average_fl_bot,3),round(Ur,3),format(raeff,'.3e'),round(dlid,3),round(Tlid,3))
print(model,Ra0/1e5,f,Ea/1e5,)


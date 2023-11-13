#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 11:48:34 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import time
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import function_savedata as fs


figpath = '/Users/chingchen/Desktop/figure/StagYY/'
labelsize = 20
bwith = 3


path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/'
modelpath = '/Users/chingchen/Desktop/model/'
model = 'Ra3.2e5_Ea1e5_f0.8'


tw2=0.17
tw1=tw2-0.06


#Ra1e4_Ea1e5_f0.75
#Ra1e4_Ea1e6
#Ra1e5_Ea1e6
#Ra1e5_Ea1e7
#Ra3.2e4_Ea1e5
#Ra3.2e4_Ea1e6_f0.7
#Ra3.2e4_Ea1e7_f0.8
#Ra3.2e5_Ea1e5_f0.8


model_information = pd.read_csv(workpath+'model_information_all.csv',sep=',')

fig_ftime = 1
choosen_time_window = 1
plot_average = 1
plot_Tprofile = 0
save= 1
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
H = model_information.H[jj]

ff = pd.read_csv(path+model+'_scaling_time_data.csv')
t0 = time.time()
print("----t0----: ", t0- start)
kk = np.arange(0.03,2.5,0.005) # time shift = 0.005
tt1=np.zeros(len(kk))
tt2=np.zeros(len(kk))
anb=np.zeros(len(kk))
ant=np.zeros(len(kk))
afb=np.zeros(len(kk))
aft=np.zeros(len(kk))

if fig_ftime:
    fig,(ax) = plt.subplots(1,1,figsize=(12,10))
    #ax.plot(ff.time,ff.Nu_bot,color='#AE6378',label = 'Nu_bot')
    #ax.plot(ff.time,ff.Nu_top,color='#35838D',label = 'Nu_top',lw=3)
    Ftop_pred = f**2 * ff.F_bot
    ax.plot(ff.time,Ftop_pred,color='#AE6378',label = 'F_bot')
    ax.plot(ff.time,ff.F_top,color='#35838D',label = 'F_top',lw=3)
    
if choosen_time_window:
    for ii,tt in enumerate(kk): 
        time_window1=tt
        time_window2=tt+0.05 # time window = 0.1
        if len(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])==0:
            break
        average_Nu_bot = np.median(ff.Nu_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        average_Nu_top = np.median(ff.Nu_top[(ff.time>time_window1) & (ff.time<time_window2)])
        average_fl_bot = np.median(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        average_fl_top = np.median(ff.F_top[(ff.time>time_window1) & (ff.time<time_window2)])
        
        tt1[ii] = np.round(time_window1,2)
        tt2[ii] = np.round(time_window1,2)
        
       
        afb[ii] = average_fl_bot
        aft[ii] = average_fl_top
    
    tt1 = tt1[afb>0]
    tt2 = tt2[afb>0]
    afb = afb[afb>0]
    if len(tt1[(afb>np.mean(afb)+np.std(afb))])==0:
        print(model)
    
    #time_window1=np.max(tt1[(anb>np.mean(anb)+np.std(anb))])+0.3 # minimum time + 0.3
    #time_window2=np.max(tt2[(anb<np.mean(anb)+np.std(anb))])-0.05 # maximum time - 0.05
    time_window1=tw1
    time_window2=tw2
    

    if plot_average:
        ax.legend(fontsize =labelsize)
        ax.set_title(model,fontsize=labelsize)
        ax.set_ylim(2,2.1)
        for aa in [ax]:
            aa.set_xlim(0,time_window2+0.05)
            aa.tick_params(labelsize=labelsize)
            aa.set_ylabel('Nu$_{bot}$',fontsize = 20)
            aa.axvline(x=time_window1,color = '#97795D')
            aa.axvline(x=time_window2,color = '#6B0D47')
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
                

# calculate the average Nu and heat flow in given time window    
fbot = np.average(ff.F_bot[(ff.time>time_window1)&(ff.time<time_window2)])
ftop = np.average(ff.F_top[(ff.time>time_window1)&(ff.time<time_window2)])

#-----------Find Average Temperature and the Stagnant Lid Thickness------------
end = int((time_window2-0.1)*1000)
new_thickness = np.zeros(end)
avg_temp = np.zeros(end)

for i in range(1,end+1):
    ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = rprof_header_list)    
    x = np.array(ff.vzabs*ff.Tmean)
    y = np.array(ff.r)
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
        avg_temp[i-1] = np.average(ff.Tmean[y<y[x==inf_point]])



rprof_ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
              sep = '\\s+',header = None,names = rprof_header_list) 
Tmm = np.average(rprof_ff.Tmean[(rprof_ff.r>0.3)*(rprof_ff.r<0.4)])
kmm = np.average(rprof_ff.tcondmean[(rprof_ff.r>0.3)*(rprof_ff.r<0.4)])
average_t = Tmm
dlid = 1-np.average(new_thickness[500:])

gamma = np.log(np.array(Ea))
rsurf = Ra0/np.exp(gamma/2)
raeff = rsurf*np.exp(gamma*Tmm) 

print('##############################################################')
print(model,round(rsurf,3),round(gamma,3),f,round(Tmm,3),round(ftop,3),
      round(fbot,3),format(raeff,'.3e'),round(time_window1,3),
      round(time_window2,3),format(Ea,'.1e'),)
print(model,Ra0/1e5,f,Ea/1e5,)
print(round(ftop,3),round(f**2 * fbot,3),round(f**2 * fbot -ftop,3))

ax.set_ylim(ftop-0.2,ftop+0.2)
#ax.set_ylim(1,3.5)
#ax2.set_ylim(ftop-0.2,ftop+0.2)
if plot_Tprofile:
    fig3,(ax6) = plt.subplots(1,1,figsize=(6,10))
    ax6.plot(rprof_ff.Tmean,rprof_ff.r,color ='#414F67',lw=3)
    ax6.axvline(x=Tmm,color ='#97795D')
    for aa in [ax6]:
        aa.set_xlim(0.85,1.05)
        aa.tick_params(labelsize=labelsize)
        aa.set_ylim(0,1)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
Ur = 0
Tlid = 0
if save:
    savearray = np.array([round(rsurf,5),f,round(Ea,0),H,round(Tmm,3),round(ftop,3),
          round(fbot,3),round(Ur,3),round(raeff,3),round(dlid,3),round(Tlid,3)])
    fs.save_1txt(model+'_scaling_grep_data',path,savearray)
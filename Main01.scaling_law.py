#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  29 16:44:37 2023

@author: chingchen
"""
import sys
import pandas as pd
import numpy as np
import time
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import function_savedata as fs


labelsize = 20
bwith = 3
local = 1

if local:
    figpath = '/Users/chingchen/Desktop/figure/StagYY/'
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    modelpath = '/Users/chingchen/Desktop/model/'
else:
    figpath = '/lfs/jiching/figure/'
    path = '/lfs/jiching/ScalingLaw_model/data_scaling/'
    workpath = '/lfs/jiching/ScalingLaw_model/'
    modelpath = '/lfs/jiching/ScalingLaw_model/23agu/'
model = sys.argv[1]

tw1 = 0.2
tw2 = 0.24
xmin,xmax=2.6,9.3
model_information = pd.read_csv(workpath+'model_information_23end.csv',sep=',')  

fig_ftime = 1
choosen_time_window = 1
plot_average = 1
plot_Tprofile = 0
save = 1


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

f = model_information.f[jj]
Ea = model_information.Ea[jj]
Ra0 = model_information.Ra0[jj]
H = model_information.H[jj]

ff = pd.read_csv(path+model+'_scaling_time_data.csv')
kk = np.arange(0.03,2.5,0.005) # time shift = 0.005
tt1=np.zeros(len(kk))
tt2=np.zeros(len(kk))
anb=np.zeros(len(kk))
ant=np.zeros(len(kk))
afb=np.zeros(len(kk))
aft=np.zeros(len(kk))


if fig_ftime:
    fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
    Ftop_pred = f**2 * ff.F_bot + (1+f+f**2)/3 *H
    ax.plot(ff.time,Ftop_pred,color='#AE6378',label = 'F_bot')
    ax.plot(ff.time,ff.F_top,color='#35838D',label = 'F_top',lw=3)
    ax.set_title(model,fontsize=labelsize)
    ax.set_ylabel('F$_{top}$',fontsize = labelsize)
    
if choosen_time_window:
    for ii,tt in enumerate(kk): 
        time_window1=tt
        time_window2=tt+0.005 # time window = 0.1
        if len(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])==0:
            break
        average_fl_top = np.median(ff.F_top[(ff.time>time_window1) & (ff.time<time_window2)])
        tt1[ii] = np.round(time_window1,2)
        tt2[ii] = np.round(time_window1,2)
        aft[ii] = average_fl_top
    
    tt1 = tt1[aft>0]
    tt2 = tt2[aft>0]
    aft = aft[aft>0]
    if len(tt1[(aft>np.mean(aft)+np.std(aft))])==0:
        print(model)
    
    time_window1=np.max(tt1[(aft>np.mean(aft)+0.1*np.std(ff.Nu_bot))])+0.001 # minimum time + 0.3
    time_window2=np.max(tt2[(aft<np.mean(aft)+0.1*np.std(ff.Nu_top))])-0.001 # maximum time - 0.05
    #time_window1 = tw1
    #time_window2 = tw2
    
    if plot_average:
        ax2.scatter(tt1,aft,color='#35838D')
        ax2.set_ylabel('F$_{bot}$',fontsize = 20)
        ax.legend(fontsize =labelsize)
        for aa in [ax,ax2]:
            aa.set_xlim(0,time_window2+0.05)
            aa.tick_params(labelsize=labelsize)
            aa.set_ylim(xmin,xmax)
            aa.axvline(x=time_window1,color = '#97795D')
            aa.axvline(x=time_window2,color = '#6B0D47')
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
            
# calculate the average Nu and heat flow in given time window    
average_fl_bot = np.average(ff.F_bot[(ff.time>time_window1)&(ff.time<time_window2)])
average_fl_top = np.average(ff.F_top[(ff.time>time_window1)&(ff.time<time_window2)])


#-----------Find Average Temperature and the Stagnant Lid Thickness------------
end = int((time_window2)*1000)
new_thickness = np.zeros(end)
avg_temp = np.zeros(end)

for i in range(200,end+1):
    rof = pd.read_csv(workpath+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = rprof_header_list)    
    x = np.array(rof.vzabs*rof.Tmean)
    y = np.array(rof.r)
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
        avg_temp[i-1] = np.average(rof.Tmean[y<y[x==inf_point]])

rprof_ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
              sep = '\\s+',header = None,names = rprof_header_list) 
Tmm = np.average(rprof_ff.Tmean[(rprof_ff.r>0.3)*(rprof_ff.r<0.5)])
average_t = Tmm
dlid = 1-np.average(new_thickness[20:])
ftop = average_fl_top
fbot = average_fl_bot
gamma = np.log(np.array(Ea))
rsurf = Ra0/np.exp(gamma/2)
raeff = rsurf*np.exp(gamma*average_t)
flid = 1-(1-f)*dlid
Tlid = dlid/flid*(ftop-1/6*H*(2-flid-flid**2)/(1-f))
Ur = (1+f+f**2)/3*H/ftop
#print(model,rsurf,gamma,average_t,1e-5,average_Nu_bot,1e-3)
print('##############################################################')
print(model,round(rsurf,5),f,format(Ea,'.1e'),H,round(Tmm,3),round(ftop,3),
      round(fbot,3),round(Ur,3),format(raeff,'.3e'),round(dlid,3),round(Tlid,3))
print(round(ftop,3),round(f**2 * fbot + (1+f+f**2)/3 *H,3),round(f**2 * fbot + (1+f+f**2)/3 *H-ftop,3))


ax.set_ylim(ftop-0.2,ftop+0.2)
ax2.set_ylim(ftop-0.2,ftop+0.2)
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

if save:
    savearray = np.array([round(rsurf,5),f,round(Ea,0),H,round(Tmm,3),round(ftop,3),
          round(fbot,3),round(Ur,3),round(raeff,3),round(dlid,3),round(Tlid,3)])
    fs.save_1txt(model+'_scaling_grep_data',path,savearray)




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 03:15:37 2023

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
model = 'ca01'

model_information = pd.read_csv(workpath+'model_information.csv',sep=',')
model_information = pd.read_csv('/Users/chingchen/Desktop/StagYY_Works/model_internal_heating.csv',sep=',')  

fig_ftime = 1
choosen_time_window = 1
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

# savedata array 
leng = len(model_information)
model_array = np.zeros(leng)
time_window1_array = np.zeros(leng)
time_window2_array = np.zeros(leng)
Ra0_list = np.zeros(leng)
Ea_list = np.zeros(leng)
f_list = np.zeros(leng)

for jj in range(len(model_information)):
    if model_information.model[jj]== model:
        break 
start = time.time()
f = model_information.f[jj]
Ea = model_information.Ea[jj]
Ra0 = model_information.Ra0[jj]
H = model_information.H[jj]

file = path+model+'_time.dat'
ff = pd.read_csv(file,sep = '\\s+')    
kk = np.arange(0.1,2.5,0.01) # time shift = 0.01
tt1=np.zeros(len(kk))
tt2=np.zeros(len(kk))
anb=np.zeros(len(kk))
ant=np.zeros(len(kk))
afb=np.zeros(len(kk))
aft=np.zeros(len(kk))

if fig_ftime:
    fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
    ax.plot(ff.time,ff.Nu_bot,color='#AE6378',label = 'Nu_bot')
    ax.plot(ff.time,ff.Nu_top,color='#35838D',label = 'Nu_top',lw=3)
    Ftop_pred = f**2 * ff.F_bot + (1+f+f**2)/3 *H
    ax2.plot(ff.time,Ftop_pred,color='#AE6378',label = 'F_bot')
    ax2.plot(ff.time,ff.F_top,color='#35838D',label = 'F_top',lw=3)
    #ax2.scatter(tt1,afb,color='#35838D')
    ax.set_title(model,fontsize=labelsize)
    ax.set_ylabel('Nu$_{bot}$',fontsize = labelsize)
    ax2.set_ylabel('F$_{bot}$',fontsize = labelsize)
    for aa in [ax,ax2]:
        aa.legend(fontsize =labelsize)
        aa.tick_params(labelsize=labelsize)
        aa.set_ylim(-5,20)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if choosen_time_window:
    for ii,tt in enumerate(kk): 
        time_window1=tt
        time_window2=tt+0.001 # time window = 0.1
        if len(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])==0:
            break
        average_Nu_bot = np.median(ff.Nu_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        average_Nu_top = np.median(ff.Nu_top[(ff.time>time_window1) & (ff.time<time_window2)])
        average_fl_bot = np.median(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        average_fl_top = np.median(ff.F_top[(ff.time>time_window1) & (ff.time<time_window2)])
        
        tt1[ii] = np.round(time_window1,2)
        tt2[ii] = np.round(time_window1,2)
        anb[ii] = average_Nu_bot
        ant[ii] = average_Nu_top
        afb[ii] = average_fl_bot
        aft[ii] = average_fl_top
    
    tt1 = tt1[anb>0]
    tt2 = tt2[anb>0]
    anb = anb[anb>0]
    afb = afb[afb>0]
    if len(tt1[(anb>np.mean(anb)+np.std(anb))])==0:
        print(model)
    
    #time_window1=np.max(tt1[(anb>np.mean(anb)+0.1*np.std(ff.Nu_bot))])+0.001 # minimum time + 0.3
    #time_window2=np.max(tt2[(anb<np.mean(anb)+0.1*np.std(ff.Nu_top))])-0.01 # maximum time - 0.05
    
    time_window1 = 0.27
    time_window2 = 1.47
    if plot_average:
        fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
        #ax.scatter(ff.time,ff.Nu_bot,color='#AE6378')
        ax.set_title(model,fontsize=labelsize)
        ax.scatter(tt1,anb,color='#35838D')
        #ax2.scatter(ff.time,ff.F_bot,color='#AE6378')
        ax2.scatter(tt1,afb,color='#35838D')
        ax.set_ylabel('Nu$_{bot}$',fontsize = 20)
        ax2.set_ylabel('F$_{bot}$',fontsize = 20)
        for aa in [ax,ax2]:
            aa.tick_params(labelsize=labelsize)
            aa.set_ylim(0,5)
            aa.axvline(x=time_window1,color = '#97795D')
            aa.axvline(x=time_window2,color = '#6B0D47')
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
                
time_window1_array[jj] = time_window1
time_window2_array[jj] = time_window2
# calculate the average Nu and heat flow in given time window    
average_Nu_bot = np.average(ff.Nu_bot[(ff.time>time_window1)&(ff.time<time_window2)])
average_Nu_top = np.average(ff.Nu_top[(ff.time>time_window1)&(ff.time<time_window2)])
average_fl_bot = np.average(ff.F_bot[(ff.time>time_window1)&(ff.time<time_window2)])
average_fl_top = np.average(ff.F_top[(ff.time>time_window1)&(ff.time<time_window2)])

#print('average_Nu_bot = ', np.round(average_Nu_bot,3))

print("total time taken this loop: ", time.time() - start)
### Find Average Temperature and the Stagnant Lid Thickness
end = int((time_window2-0.1)*1000)
lid_thickness = np.zeros(end)
new_thickness = np.zeros(end)
avg_temp = np.zeros(end)
time = np.linspace(1,end,end)
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
dlid = 1-np.average(new_thickness[200:])
Ftop = average_fl_top
Fbot = average_fl_bot
gamma = np.log(np.array(Ea))
rsurf = Ra0/np.exp(gamma/2)
raeff = rsurf*np.exp(gamma*average_t)
flid = 1-(1-f)*dlid
Tlid = dlid/flid*(Ftop-1/6*H*(2-flid-flid**2)/(1-f))
Ur = (1+f+f**2)/3*H/Ftop
#print(model,rsurf,gamma,average_t,1e-5,average_Nu_bot,1e-3)
print('##############################################################')
print(model,round(rsurf,3),f,Ea,H,round(Tmm,3),round(Ftop,3),
      round(Fbot,3),round(Ur,3),round(raeff,3),round(dlid,3),round(Tlid,3))

if save_data:
    print(model,rsurf,f,Ea,H,Tmm,Ftop,Fbot,Ur,raeff,dlid,Tlid)
    qqq=DataFrame(model_information.model,columns=['model'])
    
    n6 = pd.concat([qqq,Ra,f,Ea,aaa,bbb],axis=1)
    n6.to_csv(workpath+'model_list'+'.csv',index=False)
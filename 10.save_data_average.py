#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 01:12:20 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
plt.rcParams["font.family"] = "Times New Roman"


savepath='/Users/chingchen/Desktop/data/'
labelsize = 20
bwith = 3
save_data = 1
fig_T    = 0
fig_Nu_t = 0

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
path = '/Users/chingchen/Desktop/data/'
rprof_path='/Users/chingchen/Desktop/model/'
model_list = ['w0201','w0202','w0204','w0205','w0207','w0208',
              'w0210','w0213','w0214','w0216',
             'w0219','w0220','w0222',
              'w0801','w0802','w0803','w0804','w0805','w0806','w0807','w0808',]
#              'w0809','w0810']
#

model_list=['w0213','w0214','w0216','w0217','w0219','w0220'
,'w0203','w0206','w0209','w0211','w0215','w0212','w0218','w0221','w0223',
'w0217','w0224','w0225','w0226','w0227','w0228','w0229','w0230',
'w0811','w0812','w0813','w0814','w0815','w0816',]

model_information = pd.read_csv(path+'model_list.csv') 
# savedata array 
leng = len(model_information)
model_array = np.zeros(leng)
time_window1_array = np.zeros(leng)
time_window2_array = np.zeros(leng)
average_Nu_bot_array = np.zeros(leng)
average_Nu_top_array = np.zeros(leng)
average_fl_bot_array = np.zeros(leng)
average_fl_top_array = np.zeros(leng)
average_t_array = np.zeros(leng)
lid_list = np.zeros(leng)
Ra0_list = np.zeros(leng)
Ea_list = np.zeros(leng)
f_list = np.zeros(leng)

for jj in range(len(model_information)):
    model = model_information.model[jj]
    f = model_information.f[jj]
    Ea = model_information.Ea[jj]
    Ra0 = model_information.Ra0[jj]
    time_window1 = model_information.time_window1[jj]
    time_window2 = model_information.time_window2[jj]
    
    file = path+model+'_time.dat'
    ff = pd.read_csv(file,sep = '\\s+')    
    average_Nu_bot = np.average(ff.Nu_bot[(ff.time>time_window1)*(ff.time<time_window2)])
    average_Nu_top = np.average(ff.Nu_top[(ff.time>time_window1)*(ff.time<time_window2)])
    average_fl_bot = np.average(ff.F_bot[(ff.time>time_window1)*(ff.time<time_window2)])
    average_fl_top = np.average(ff.F_top[(ff.time>time_window1)*(ff.time<time_window2)])
    Tmm = np.zeros(abs(int(time_window1*100)-int(time_window2*100)))
    new_thickness = np.zeros(abs(int(time_window1*100)-int(time_window2*100)))
    for ss,time in enumerate(range(int(time_window1*100),int(time_window2*100))):
        rprof_ff = pd.read_csv(rprof_path+model+'/datafile/'+model+'_data_'+str(time)+'.txt',
                      sep = '\\s+',header = None,names = rprof_header_list) 
        Tmm[ss] = np.average(rprof_ff.Tmean[(rprof_ff.r>0.4)*(rprof_ff.r<0.5)])
        x = np.array(rprof_ff.vzabs*rprof_ff.Tmean)
        y = np.array(rprof_ff.r)
        smooth_d2 = np.gradient(np.gradient(x))
        infls = np.where(np.diff(np.sign(smooth_d2)))[0]
        if len(x[infls])!=0:
            inf_point = x[infls][-1]
            index_inf_point = [i for i, kk in enumerate(x) if kk == inf_point][0]
            a = (y[index_inf_point]-y[index_inf_point+1])/(x[index_inf_point]-x[index_inf_point+1])
            b = y[infls][-1]-a*inf_point
            line_x = np.linspace(0,2400)
            line_y = a*line_x + b
            new_thickness[ss] = line_y[0]
    average_t = np.average(Tmm)   
    lid = np.average(new_thickness)

    if fig_T :
        fig3,(ax3) = plt.subplots(1,1,figsize=(12,6))
        kk = 1
        ax3.set_title(model, fontsize = labelsize)
        ax3.plot(ff.time, ff.Tmean,color = newcolors[kk],lw=5,label = 'Tmean')
        # ax3.plot(ff.time, ff.Tmean,color = newcolors[kk+1],lw=3,label = 'F_bot')
        
        ax3.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
            ax3.spines[axis].set_linewidth(bwith)
        ax3.grid()
        ax3.set_xlim(0,2.5)
        ax3.set_ylim(0,1)
        ax3.set_xlabel('time',fontsize = labelsize)
        ax3.set_ylabel('Temperature',fontsize = labelsize)
        ax3.legend(fontsize = 20)
        ax3.set_title(model,fontsize = labelsize)
        ax3.axvspan(time_window1,time_window2, facecolor='#2ca02c')
        ax3.axvspan(time_window1,time_window2, facecolor='#2ca02c')
    if fig_Nu_t:
        fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
        kk = 1
        ax.plot(ff.time, ff.Nu_bot,color = newcolors[kk+1],lw=3,label = 'Nu_bot')   
        ax.plot(ff.time, ff.Nu_top,color = newcolors[kk],lw=3,label = 'Nu_top')
            
        ax.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(bwith)
        ax.grid()
        ax.set_xlim(0,2.5)
        ax.set_ylim(2,8)
        ax.set_ylabel('heat flux',fontsize = labelsize)
        ax.legend(fontsize = 20)
        ax2.plot(ff.time, ff.F_bot*(f**2),color = newcolors[kk+1],lw=3,label = 'F_bot')
        ax2.plot(ff.time, ff.F_top,color = newcolors[kk],lw=3,label = 'F_top')
        ax2.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
            ax2.spines[axis].set_linewidth(bwith)
        ax2.grid()
        ax2.set_xlim(0,2.5)
        ax2.set_ylim(2,5)
        ax2.set_xlabel('time',fontsize = labelsize)
        ax2.set_ylabel('heat flux',fontsize = labelsize)
        ax2.legend(fontsize = 20)
        ax.set_title(model,fontsize = labelsize)
        ax.axvspan(time_window1,time_window2, facecolor='#2ca02c')
        ax2.axvspan(time_window1,time_window2, facecolor='#2ca02c')
    if save_data:
     
        average_Nu_bot_array[jj] = average_Nu_bot 
        average_Nu_top_array[jj] = average_Nu_top 
        average_fl_bot_array[jj] = average_fl_bot
        average_fl_top_array[jj] = average_fl_top
        average_t_array[jj]=average_t
        lid_list[jj] = lid


if save_data:         
    moodel=DataFrame(model_information.model,columns=['model'])
    f=DataFrame(model_information.f,columns=['f'])
    Ea=DataFrame(model_information.Ea,columns=['Ea'])
    Ra=DataFrame(model_information.Ra0,columns=['Ra0'])
    tt1=DataFrame(model_information.time_window1,columns=['time_window1'])
    tt2=DataFrame(model_information.time_window2,columns=['time_window2'])
    lid=DataFrame(lid_list,columns=['lid'])
    aaa=DataFrame(average_Nu_bot_array,columns=['average_Nu_bot'])
    bbb=DataFrame(average_Nu_top_array,columns=['average_Nu_top'])
    ccc=DataFrame(average_fl_bot_array,columns=['average_fl_bot'])
    ddd=DataFrame(average_fl_top_array,columns=['average_fl_top'])
    fff=DataFrame(average_t_array,columns=['Tm'])
    n6 = pd.concat([moodel,f,Ra,Ea,tt1,tt2,lid,aaa,bbb,ccc,ddd,fff],axis=1)
    n6.to_csv(path+'average_data'+'.csv',index=False)
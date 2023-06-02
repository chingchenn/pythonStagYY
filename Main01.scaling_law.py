#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 14:22:28 2023

@author: chingchen
"""

import os
import numpy as np
import pandas as pd
import function_savedata as fs
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame


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

### SETTING ###

time_shift = 0.01
time_window = 0.1

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
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

plt.rcParams["font.family"] = "Times New Roman"
model_information = pd.read_csv(workpath+'model_information.csv',sep=',') 
labelsize = 20
bwith = 3
scatter_size = 100


### DO WHAT ###

split_dat = 1
save_data = 1
plot_average = 0
fig_compare_T = 0
fig_compare_Fu = 0
fig_compare_Nu = 0
running_average = 1
plot_scaling_result = 0


if split_dat:
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        file = path+model+'_rprof.dat'
        sourceFileName = file
        sourceFileData = open(file,'r')
        ListOfLine = sourceFileData.read().splitlines()
        new_data = ListOfLine[1:]
        n = len(new_data)
        p = 129
        m = int(n/p)
        for i in range(m):
            destFileName = modelpath+model+'/datafile/'
            if not os.path.isdir(modelpath+model):
                os.mkdir(modelpath+model)
            if not os.path.isdir(destFileName):
                os.mkdir(destFileName)
            destFileData = open(destFileName+model+'_data_'+str(i)+'.txt','w')
            
            if (i==m-1):
                for line in new_data[i*p+1:]:
                    destFileData.write(line+'\n')
            else:
                for line in new_data[i*p+1:(i+1)*p]:
                    destFileData.write(line+'\n')
            destFileData.close()
        print(model+'==DONE==',m-1)
if running_average: # Using the moving window to choose the time window 
    # savedata array (model number length)
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
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        f = model_information.f[jj]
        Ea = model_information.Ea[jj]
        Ra0 = model_information.Ra0[jj]
        file = path+model+'_time.dat'
        ff = pd.read_csv(file,sep = '\\s+') 
        # time series array (window length)
        window = np.arange(0.0,2.5,time_shift) # time_shift = 0.01
        tt1 = np.zeros(len(window))
        tt2 = np.zeros(len(window))
        anb = np.zeros(len(window))
        afb = np.zeros(len(window))
        for ii,tt in enumerate(window): 
            time_window1 = tt
            time_window2 = tt + time_window # time_window = 0.1
            if len(ff.F_bot[(ff.time>time_window1) * (ff.time<time_window2)])==0:
                break
            average_Nu_bot = np.median(ff.Nu_bot[(ff.time>time_window1) * (ff.time<time_window2)])
            average_fl_bot = np.median(ff.F_bot[(ff.time>time_window1) * (ff.time<time_window2)])
            tt1[ii] = np.round(time_window1,2)
            tt2[ii] = np.round(time_window1,2)
            anb[ii] = average_Nu_bot
            afb[ii] = average_fl_bot
    
        tt1 = np.array(tt1)
        tt2 = np.array(tt2)
        if len(tt1[(anb>np.mean(anb)+np.std(anb))])==0:
            print('========== '+model+' Cannot find suitble time window ==========')
            continue
        time_window1=np.max(tt1[(anb>np.mean(anb)+np.std(anb))])+0.01 # minimum time + 0.01
        time_window2=np.max(tt2[(anb<np.mean(anb)+np.std(anb))]) # maximum time - 0.5
        if plot_average:
            fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
            ax.set_title(model,fontsize=labelsize)
            ax.scatter(tt1,anb,color='#35838D')
            ax2.scatter(tt1,afb,color='#35838D')
            for aa in [ax,ax2]:
                aa.tick_params(labelsize=labelsize)
                aa.set_xlim(0,2.5)
                for axis in ['top','bottom','left','right']:
                    aa.spines[axis].set_linewidth(bwith)
    
        
#       # calculate the average Nu and heat flow in given time window    
        average_Nu_bot = np.average(ff.Nu_bot[(ff.time>time_window1)*(ff.time<time_window2)])
        average_Nu_top = np.average(ff.Nu_top[(ff.time>time_window1)*(ff.time<time_window2)])
        average_fl_bot = np.average(ff.F_bot[(ff.time>time_window1)*(ff.time<time_window2)])
        average_fl_top = np.average(ff.F_top[(ff.time>time_window1)*(ff.time<time_window2)])
        ### Find Average Temperature and the Stagnant Lid Thickness
        Tmm = np.zeros(abs(int(time_window1*100)-int(time_window2*100)))
        new_thickness = np.zeros(abs(int(time_window1*100)-int(time_window2*100)))
        for ss,time in enumerate(range(int(time_window1*100),int(time_window2*100))):
            rprof_ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(time)+'.txt',
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
    
        if save_data:
            time_window1_array[jj] = time_window1
            time_window2_array[jj] = time_window2
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
        tt1=DataFrame(time_window1_array,columns=['time_window1'])
        tt2=DataFrame(time_window2_array,columns=['time_window2'])
        lid=DataFrame(lid_list,columns=['lid'])
        aaa=DataFrame(average_Nu_bot_array,columns=['average_Nu_bot'])
        bbb=DataFrame(average_Nu_top_array,columns=['average_Nu_top'])
        ccc=DataFrame(average_fl_bot_array,columns=['average_fl_bot'])
        ddd=DataFrame(average_fl_top_array,columns=['average_fl_top'])
        fff=DataFrame(average_t_array,columns=['Tm'])
        ra = np.array(model_information.Ra0)
        gamma = np.log(np.array(model_information.Ea))
        nu = average_Nu_bot_array
        xx1 = np.ones(len(nu))*1e-5
        xx2 = np.ones(len(nu))*1e-3
        thetam = (average_t_array)
        rsurf = ra/np.exp(gamma/2)
        fs.save_6txt('scalinglaw_inversion',path,rsurf,gamma,thetam,xx1,nu,xx2)
#        n6 = pd.concat([moodel,f,Ra,Ea,tt1,tt2,lid,aaa,bbb,ccc,ddd,fff],axis=1)
#        n6.to_csv(path+'average_data'+'.csv',index=False)
x=[]
y=[]
if plot_scaling_result:
    fig5,(ax7) = plt.subplots(1,1,figsize=(7,7))
    # fig8,(ax8) = plt.subplots(1,1,figsize=(7,7))
    model_data = pd.read_csv(path+'average_data'+'.csv') 
    for jj in range(len(model_data)):
        model = model_data.model[jj]
        f = model_data.f[jj]
        Ea = model_data.Ea[jj]
        Ra0 = model_data.Ra0[jj]
        time_window1 = model_data.time_window1[jj]
        time_window2 = model_data.time_window2[jj]
        average_Nu_bot = model_data.average_Nu_bot[jj]
        average_fl_bot = model_data.average_fl_bot[jj]
        average_t = model_data.Tm[jj]
        
        ## Setting shape of points to distingish the f 
        # if f==0.3:
        #     marker ="^"
        # elif f==0.5:
        #     marker ="s"
        if f==0.3 or f == 0.5 or f ==0.7 or f ==0:
            continue
        # elif f==0.7:
        #     marker='o'
        elif f == 0.8:
            marker='p'
        ## Setting colour of points to distingish the Ea
        if Ea ==1e5:
            kk = 2
            label='1e5'
        elif Ea ==1e6:
            kk = 0
            label='1e6'
        elif Ea ==1e7:
            kk = 8
            label='1e7'
        elif Ea ==1e8:
            kk = 1
            label='1e8'    

        plt.rcParams['font.size'] = 20
        if average_Nu_bot>0:
            print(model,Ra0,Ea/1e5,average_Nu_bot)
            if fig_compare_Fu :
                if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210':
                    ax7.scatter(Ra0,average_fl_bot*(f**2),color=newcolors[kk+1],label = label,s=scatter_size,marker=marker)
                else:
                    ax7.scatter(Ra0,average_fl_bot*(f**2),color=newcolors[kk+1],s=scatter_size,marker=marker)
                ax7.set_ylabel('F_bot',fontsize = labelsize)
                ax7.set_yscale('log')
                
            if fig_compare_Nu :
                if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210' or model=='w0812' or model=='w0811' or model=='w0810' or model=='w0809':
                    ax7.scatter(Ra0,average_Nu_bot,color=newcolors[kk+1],label = label,s=scatter_size,marker=marker)
                else:
                    ax7.scatter(Ra0,average_Nu_bot,color=newcolors[kk+1],s=scatter_size,marker=marker)
                ax7.set_ylabel('Nu',fontsize = labelsize)
                ax7.set_yscale('log')
                
            if fig_compare_T :
                if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210':
                    ax7.scatter(Ra0,average_t,color=newcolors[kk+1],label = label,s=scatter_size,marker=marker)
                else:
                    ax7.scatter(Ra0,average_t,color=newcolors[kk+1],s=scatter_size,marker=marker)
                ax7.set_ylabel('Temp',fontsize = labelsize)
        
            x.append(Ra0)
            y.append(average_Nu_bot)
        # if Ra0==1e5:
        #     kk = 2
        #     label='1e5'
        # elif Ra0==3.2e4:
        #     kk = 0
        #     label='3.2e4'
        # elif Ra0==5.6e4:
        #     kk = 10
        #     label='5.6e4'
        # elif Ra0 ==1e4:
        #     kk = 4
        #     label='1e4' 
        # else:
        #     kk = 8
        #     label = 'nuknow'
        # if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210':
        #     ax8.scatter(Ea,average_t,color=newcolors[kk+1],label = label,s=scatter_size,marker=marker)
        # else:
        #     ax8.scatter(Ea,average_t,color=newcolors[kk+1],s=scatter_size,marker=marker)

        for axis in ['top','bottom','left','right']:
            ax7.spines[axis].set_linewidth(bwith)
            # ax8.spines[axis].set_linewidth(bwith)
    
        # ax7.grid()
        ax7.set_xscale('log')
        # ax7.set_yscale('log')
        
        ax7.legend(fontsize=labelsize)
        ax7.set_xlabel('Ra0',fontsize = labelsize)
        
        # ax7.set_xticks(np.logspace(3,5,3))
        # ax7.set_yticks(np.arange(0, 20, 1))
        ax7.tick_params(axis='both', which='major',labelsize=labelsize)
        # ax7.set_ylim(1,10)

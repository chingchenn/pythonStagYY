#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 12:16:33 2023

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
plot_average = 1
save_data = 0
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/'
model = 's01'


model_information = pd.read_csv(workpath+'model_information.csv',sep=',')
model_information = pd.read_csv('/Users/chingchen/Desktop/StagYY_Works/model_internal_heating_tem.csv',sep=',')  

# savedata array 
leng = len(model_information)
model_array = np.zeros(leng)
time_window1_array = np.zeros(leng)
time_window2_array = np.zeros(leng)
Ra0_list = np.zeros(leng)
Ea_list = np.zeros(leng)
f_list = np.zeros(leng)
for jj in range(len(model_information)):
    start = time.time()
    model = model_information.model[jj]
    f = model_information.f[jj]
    Ea = model_information.Ea[jj]
    Ra0 = model_information.Ra0[jj]
    
    file = path+model+'_time.dat'
    ff = pd.read_csv(file,sep = '\\s+')    
    kk = np.arange(0.2,2.5,0.01) # time shift = 0.02
    tt1=np.zeros(len(kk))
    tt2=np.zeros(len(kk))
    anb=np.zeros(len(kk))
    ant=np.zeros(len(kk))
    afb=np.zeros(len(kk))
    aft=np.zeros(len(kk))
    
    # tt1=[]
    # tt2=[]
    # anb=[]
    # ant=[]
    # afb=[]
    # aft=[]
    for ii,tt in enumerate(kk): 
        time_window1=tt
        time_window2=tt+0.1 # time window = 0.1
        if len(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])==0:
            break
        average_Nu_bot = np.median(ff.Nu_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        #average_Nu_top = np.median(ff.Nu_top[(ff.time>time_window1) & (ff.time<time_window2)])
        average_fl_bot = np.median(ff.F_bot[(ff.time>time_window1) & (ff.time<time_window2)])
        #average_fl_top = np.median(ff.F_top[(ff.time>time_window1) & (ff.time<time_window2)])
        
        tt1[ii] = np.round(time_window1,2)
        tt2[ii] = np.round(time_window1,2)
        anb[ii] = average_Nu_bot
        #ant[ii]=average_Nu_top
        afb[ii] = average_fl_bot
        #aft[ii] = average_fl_top
        
        # tt1.append(np.round(time_window1,2))
        # tt2.append(np.round(time_window1,2))
        # anb.append(average_Nu_bot)
        # #ant.append(average_Nu_top)
        # afb.append(average_fl_bot)
        # #aft.append(average_fl_top)
    
    tt1 = tt1[anb>0]
    tt2 = tt2[anb>0]
    anb = anb[anb>0]
    afb = afb[afb>0]
    #timemid = time.time()
    #print("total time taken this loop: ", timemid - start)
    if len(tt1[(anb>np.mean(anb)+np.std(anb))])==0:
        print(model)
        continue
    time_window1=np.max(tt1[(anb>np.mean(anb)+np.std(anb))])+0.3 # minimum time + 0.3
    time_window2=np.max(tt2[(anb<np.mean(anb)+np.std(anb))])-0.05 # maximum time - 0.05
    print(model,f,Ea,Ra0,time_window1,time_window2)
    if plot_average:
        fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,10))
        #ax.scatter(ff.time,ff.Nu_bot,color='#AE6378')
        ax.set_title(model,fontsize=labelsize)
        ax.scatter(tt1,anb,color='#35838D')
        #ax2.scatter(ff.time,ff.F_bot,color='#AE6378')
        ax2.scatter(tt1,afb,color='#35838D')

        for aa in [ax,ax2]:
            aa.tick_params(labelsize=labelsize)
            #aa.set_ylim(0,10)
            aa.axvline(x=time_window1,color = 'b')
            aa.axvline(x=time_window2,color = 'r')
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
    # if save_data:
    time_window1_array[jj] = time_window1
    time_window2_array[jj] = time_window2
    Ra0_list[jj] = Ra0
    Ea_list[jj] = Ea
    f_list[jj] = f
    print("total time taken this loop: ", time.time() - start)
if save_data:
    qqq=DataFrame(model_information.model,columns=['model'])
    f=DataFrame(f_list,columns=['f'])
    Ea=DataFrame(Ea_list,columns=['Ea'])
    Ra=DataFrame(Ra0_list,columns=['Ra0'])
    aaa=DataFrame(time_window1_array,columns=['time_window1'])
    bbb=DataFrame(time_window2_array,columns=['time_window2'])
    n6 = pd.concat([qqq,Ra,f,Ea,aaa,bbb],axis=1)
    n6.to_csv(workpath+'model_list'+'.csv',index=False)
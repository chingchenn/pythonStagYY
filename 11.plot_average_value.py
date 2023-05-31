#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 23:51:41 2023

@author: chingchen
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
figpath = '/Users/chingchen/Desktop/figure/StagYY/'
path = '/Users/chingchen/Desktop/data/'
labelsize = 12
bwith = 3
fig_T    = 0
fig_mobility = 0
fig_compare_Fu = 0
fig_compare_Nu = 1
fig_compare_T = 0
scatter_size = 100

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']

model_list = [
'w0201','w0204','w0207','w0210',
'w0801','w0802','w0803','w0804',
'w0213','w0216','w0219','w0222',
'w0225','w0228',] 

model_list = ['w0201','w0202','w0204','w0205',
              'w0207','w0208','w0210','w0213',
              'w0214','w0216','w0219','w0220',
              'w0222','w0801','w0802','w0803',
              'w0804','w0805','w0806','w0807','w0808',]


#'w0805','w0806','w0807','w0808',]
               # 'w0813','w0814','w0815','w0816',]
# model_list = ['w0201','w0202','w0204','w0205',
#               'w0207','w0208','w0210','w0213',
#               'w0214','w0216','w0219','w0220',
#               'w0222','w0801','w0802','w0803',
#               'w0804','w0805','w0806','w0807',
#               'w0808','w0809','w0810','w0225','w0228',
# ]
fig5,(ax7) = plt.subplots(1,1,figsize=(12,12))
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
    ## Setting shape of points to distingish the f 
    # if f==0.3:
    #     marker ="^"
    # elif f==0.5:
    #     marker ="s"
    if f==0.3 or f == 0.5:
        continue
    elif f==0.7:
        marker='o'
    elif f == 0.8:
        marker='p'
    
    if fig_compare_Fu :
        if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210':
            ax7.scatter(Ra0,average_fl_bot*(f**2),color=newcolors[kk+1],label = label,s=scatter_size,marker=marker)
        else:
            ax7.scatter(Ra0,average_fl_bot*(f**2),color=newcolors[kk+1],s=scatter_size,marker=marker)
        ax7.set_ylabel('F_bot',fontsize = labelsize)
        ax7.set_yscale('log')
    if fig_compare_Nu :
        if model=='w0201' or model=='w0204' or model=='w0207'or model=='w0210':
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

for axis in ['top','bottom','left','right']:
    ax7.spines[axis].set_linewidth(bwith)

# ax7.grid()
ax7.set_xscale('log')
ax7.legend(fontsize=labelsize)
ax7.set_xlabel('Ra0',fontsize = labelsize)
ax7.tick_params(labelsize=labelsize)
ax7.set_xticks(np.logspace(3,5,3))
# ax7.set_yticks(np.arange(0, 20, 1))
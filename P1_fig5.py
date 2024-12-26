#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 22:14:30 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from plot_time_vs_pow import hat_graph


labelsize = 18
bwith = 3

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#3CB371','#4682B4','#CD5C5C','#97795D','#414F67','#4198B9','#3CB371']


header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

# 
fig,(aa1,aa2) = plt.subplots(2,2,figsize=(18,14))
ax=aa1[0]
ax2=aa2[0]
axx=aa1[1]
axx3=aa2[1]
# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d13_P1.2TW_1.5wt%_D5.0km-NH3',# viscosity
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.2TW_1.5wt%_D5.0km-NH3', 
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P1.2TW_1.5wt%_D5.0km-NH3',
              ]

label_list=['1e13 Pa s','3.2e13 Pa s','5.6e13 Pa s','1e14 Pa s',]
label_list=['10$^{13}$','10$^{13.5}$','10$^{13.75}$','10$^{14}$','10$^{14.5}$']

min_zbot = []
max_zbot = []
amplitude_zbot=[]
for i, model in enumerate(model_list):
    i = i
    # data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
    #     header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
    #------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(zbot_conv)
    mins, _  = find_peaks(zbot_conv*-1)
    if len(zbot_conv[peaks])>20:
        amplitude_zbot.append(np.median(zbot_conv[peaks][10:19]-zbot_conv.data[mins][10:19]))
    else:
        print('---- HELP ----')      
    mask_cond = data.conv[data.time_Gyr>3.5]
    mask_conv = ~ma.array(data.zbot[data.time_Gyr>3.5], mask = data.conv[data.time_Gyr>3.5]).mask
    
    zbot_cond = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_cond)
    zbot_conv = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_conv)
    min_zbot.append(np.min(zbot_conv))
    max_zbot.append(np.max(zbot_conv))
   
playerA = np.array(max_zbot)
playerB = np.array(min_zbot)
# playerA = np.array(avg_amp)
# playerB = np.zeros(len(playerA))
hat_graph(axx, label_list, [playerA, playerB], ['Player A', 'Player B'])
axx.set_ylim(161,0)
axx.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
axx.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=10)
# ---------------------------------------- figure 2 --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d13_P0.6TW_1.5wt%_D5.0km-NH3',# viscosity
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.6TW_1.5wt%_D5.0km-NH3', 
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P0.6TW_1.5wt%_D5.0km-NH3',
              ]

# model_list = ['Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.05', 
#               'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.1',
#               'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.15',
#               'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.2', 
#               'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.3',
#               'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3_core0.5',
#               ]
min_zbot = []
max_zbot = []
amplitude_zbot=[]

for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax2.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    ax2.plot(x,zbot_cond,color=colors[i],linestyle='dashed')    
    #------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(zbot_conv)
    mins, _  = find_peaks(zbot_conv*-1)
    if len(zbot_conv[peaks])>20:
        amplitude_zbot.append(np.median(zbot_conv[peaks][10:19]-zbot_conv.data[mins][10:19]))
    else:
        print('---- HELP ----')      
    mask_cond = data.conv[data.time_Gyr>3.5]
    mask_conv = ~ma.array(data.zbot[data.time_Gyr>3.5], mask = data.conv[data.time_Gyr>3.5]).mask
    
    zbot_cond = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_cond)
    zbot_conv = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_conv)
    pp1 = round(np.min(zbot_conv),1)
    pp2 = round(np.max(zbot_conv),1)
    min_zbot.append(pp1)
    max_zbot.append(pp2)
#  ------------------------------ figure setting ------------------------------
ax.set_ylim(161,0)
ax2.set_ylim(161,0)
ax2.legend(fontsize=labelsize-3)

for aa in [ax,ax2]:
    aa.set_xlabel('Time (Gyr)',fontsize=labelsize)
    aa.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
playerA = np.array(max_zbot)
playerB = np.array(min_zbot)
# playerA = np.array(avg_amp)
# playerB = np.zeros(len(playerA))
hat_graph(axx3, label_list, [playerA, playerB], ['Player A', 'Player B'])
axx3.set_ylim(161,0)
axx3.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
axx3.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=10)
for aa in [axx,axx3]:
    aa.set_xlabel('Reference viscosity (Pa s)',fontsize=labelsize)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    aa.grid()

# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figure5_v5.pdf')
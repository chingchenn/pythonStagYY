#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:49:41 2024

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


labelsize = 20
bwith = 3

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']
colors2=['#B22222','#FF8C00','#4B0082','#C71585','#D2691E']
colors=['#282130','#849DAB','#CD5C5C','#35838D',]

header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

# 
fig,(aa1,aa2) = plt.subplots(2,2,figsize=(18,12))
ax2=aa1[1]
axx4=aa2[1]
# axx3=aa3[1]

ax=aa1[0]
ax3=aa2[0]
# axx=aa3[0]
# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx5%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx20%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx50%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',# composition
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.6TW_1.5wt%_D5.0km-NH3',# power
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.8TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              ]

label_list = ['5%','20%','50%']
label_list = ['10$^{13.5}$','10$^{13.75}$','10$^{14}$']
label_list = ['0.6 TW','0.8 TW','1.0 TW']

min_zbot = []
max_zbot = []
amplitude_zbot=[]
for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.dlid, mask = mask_cond)
    zbot_conv = ma.array(data.dlid, mask = mask_conv)
    if i < 3:
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    # ax.plot(x,zbot_cond,color='b',linestyle='dashed', lw=1)
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Tm, mask = mask_cond)
    zbot_conv = ma.array(data.Tm, mask = mask_conv)
    if i < 3:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    
# ---------------------------------------- figure 2 --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx5%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx20%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx50%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.2TW_1.5wt%_D5.0km-NH3',# composition
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.8TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              ]
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
    zbot_cond = ma.array(data.dlid, mask = mask_cond)
    zbot_conv = ma.array(data.dlid, mask = mask_conv)
    
    if i < 3:
        ax2.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax2.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    # ax2.plot(x,zbot_cond,color='b',linestyle='dashed', lw=1)    
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Tm, mask = mask_cond)
    zbot_conv = ma.array(data.Tm, mask = mask_conv)
    if i < 3:
        axx4.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        axx4.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    # axx4.plot(x,zbot_cond,color='b',linestyle='dashed', lw=1)
#  ------------------------------ figure setting ------------------------------
ax.set_ylim(20,0)
ax2.set_ylim(20,0)
ax3.set_ylim(240,280)
axx4.set_ylim(240,280)
ax2.legend(fontsize=labelsize)
ax2.set_ylabel('stagnant thickness (km)',fontsize = labelsize)
ax.set_ylabel('stagnant thickness (km)',fontsize = labelsize)
ax3.set_ylabel('T$_m$ (K)',fontsize = labelsize)
axx4.set_ylabel('T$_m$ (K)',fontsize = labelsize)
# axx3.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
# axx.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
for aa in [ax,ax3,ax2,axx4]:
    aa.set_xlabel('Time (Gyr)',fontsize=labelsize)
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=5)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# for aa in [axx,axx3]:
#     for axis in ['top','bottom','left','right']:
#         aa.spines[axis].set_linewidth(bwith)
#     aa.grid()
# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figureS4_v1.pdf')
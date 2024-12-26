#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 11:58:06 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt


labelsize = 20
bwith = 3

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']
colors=['#B22222','#849DAB','#4B0082','#C71585','#D2691E']

header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

# 
fig,(aa1) = plt.subplots(1,1,figsize=(11,8))
ax=aa1

# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3_core0.05',
              'Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3_core0.1',
              'Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3_core0.15',
              'Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3_core0.2'
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx5%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx20%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx50%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx50%_eta1.0d14_P0.8TW_3.0wt%_D5.0km-NH3_Hvar',
              'Europa-tidal5_period0.14Gyr_emx50%_eta1.0d14_P0.8TW_3.0wt%_D5.0km-NH3_Hvar_2']
model_list = ['Europa-tidal5_period0.14Gyr_emx5%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx20%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx50%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx5%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3'
              ]
label_list = ['5%','10%','20%','30%','50%']
label_list = ['161 km','138 km','10%','15%','20%']

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
    zbot_cond = ma.array(data.Prad+data.Pint, mask = mask_cond)
    zbot_conv = ma.array(data.Prad+data.Pint, mask = mask_conv)
    
    ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    # ax.plot(x,zbot_cond,color='darkgreen',linestyle='dashed')
#  ------------------------------ figure setting ------------------------------
ax.set_ylim(0,2.5)
ax.set_ylabel('power',fontsize = labelsize)
for aa in [ax]:
    aa.minorticks_on()
    aa.set_xlabel('Time (Gyr)',fontsize=labelsize)
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figureS1_v1.pdf')
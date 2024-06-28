#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 15:35:44 2024

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

def hat_graph(ax, xlabels, values, group_labels):
    """
    Create a hat graph.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes to plot into.
    xlabels : list of str
        The category names to be displayed on the x-axis.
    values : (M, N) array-like
        The data values.
        Rows are the groups (len(group_labels) == M).
        Columns are the categories (len(xlabels) == N).
    group_labels : list of str
        The group labels displayed in the legend.
    """

    def label_bars(heights, rects):
        """Attach a text label on top of each bar."""
        for height, rect in zip(heights, rects):
            ax.annotate(f'{height}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 2),  # 4 points vertical offset.
                        textcoords='offset points',
                        ha='center', va='bottom',fontsize=14)

    values = np.asarray(values)
    x = np.arange(values.shape[1])
    ax.set_xticks(x, labels=xlabels)
    spacing = 0.3  # spacing between hat groups
    width = (1 - spacing) / values.shape[0]
    heights0 = values[0]
    for i, (heights, group_label) in enumerate(zip(values, group_labels)):
        style = {'fill': False} if i == 0 else {'edgecolor': 'black'}
        rects = ax.bar(x, heights - heights0,
                       width, bottom=heights0, label=group_label, **style)
        label_bars(heights, rects)
labelsize = 20
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
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.1TW_2.0wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.3TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.5TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.6TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.8TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_2.0wt%_D5.0km-NH3',]

model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.1TW_2.0wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.3TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.5TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.8TW_2.0wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_2.0wt%_D5.0km-NH3',]

model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.1TW_1.5wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.2TW_1.5wt%_D5.0km-NH3',]

model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P0.6TW_1.5wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P0.8TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.2TW_1.5wt%_D5.0km-NH3',]
label_list=['0.1 TW','0.3 TW','0.5 TW','0.6 TW','0.8 TW','1.0 TW']#,'P = 1.2 Tw','P = 1.4 Tw']
label_list=['P = 0.1 TW','P = 0.6 TW','P = 1.0 TW','P = 1.2 TW']
label_list=['0.6 TW','0.8 TW','1.0 TW','1.2 TW']#,'P = 1.4 Tw']

# rainbow = cm.get_cmap('winter',len(model_list))
# colors = rainbow(np.linspace(0, 1, len(model_list)+1))
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
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
    # fig2,(am) = plt.subplots(1,1,figsize=(8,8))
    # am.plot(x,data['melt'],color='orange',lw=2)
    # am.set_title(model)
    #------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(zbot_conv)
    mins, _  = find_peaks(zbot_conv*-1)
    if len(zbot_conv[peaks])>20:
        amplitude_zbot.append(np.median(zbot_conv[peaks][10:19]-zbot_conv.data[mins][10:19]))
    else:
        print('---- HELP ----')      
    mask_cond = data.conv[data.time_Gyr>1.5]
    mask_conv = ~ma.array(data.zbot[data.time_Gyr>1.5], mask = data.conv[data.time_Gyr>1.5]).mask
    
    zbot_cond = ma.array(data.zbot[data.time_Gyr>1.5], mask = mask_cond)
    zbot_conv = ma.array(data.zbot[data.time_Gyr>1.5], mask = mask_conv)
    min_zbot.append(np.min(zbot_conv))
    max_zbot.append(np.max(zbot_conv))
    
# fig2,(axx,axx3) = plt.subplots(2,1,figsize=(10,12))
playerA = np.array(max_zbot)
playerB = np.array(min_zbot)
hat_graph(axx, label_list, [playerA, playerB], ['Player A', 'Player B'])
axx.set_ylim(161,0)
axx.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
axx.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=10)

model_list = ['Europa-tidal1_eta3.2d13_P0.1TW_1.5wt%-NH3', # power
                'Europa-tidal1_eta3.2d13_P0.6TW_1.5wt%-NH3',
                'Europa-tidal1_eta3.2d13_P1.0TW_1.5wt%-NH3',
                'Europa-tidal1_eta3.2d13_P1.2TW_1.5wt%-NH3',]
model_list = ['Europa-tidal1_eta5.6d13_P0.6TW_1.5wt%-NH3', # power
              'Europa-tidal1_eta5.6d13_P0.8TW_1.5wt%-NH3',
              'Europa-tidal1_eta5.6d13_P1.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta5.6d13_P1.2TW_1.5wt%-NH3',]


for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax.plot(x,zbot_conv,color=colors[i],lw=3)
    # ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')  
    
# ---------------------------------------- figure 2 --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.1TW_1.5wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.1TW_1.5wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_1.5wt%_D5.0km-NH3',
              ]
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P0.6TW_1.5wt%_D5.0km-NH3', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P0.8TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P1.2TW_1.5wt%_D5.0km-NH3',
              ]

# rainbow = cm.get_cmap('winter',len(model_list))
# colors = rainbow(np.linspace(0, 1, len(model_list)+1))
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
    # fig2,(am) = plt.subplots(1,1,figsize=(8,8))
    # am.plot(x,data['melt'],color='orange',lw=2)
    # am.set_title(model)
    #------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(zbot_conv)
    mins, _  = find_peaks(zbot_conv*-1)
    if len(zbot_conv[peaks])>20:
        amplitude_zbot.append(np.median(zbot_conv[peaks][10:19]-zbot_conv.data[mins][10:19]))
    else:
        print('---- HELP ----')      
    mask_cond = data.conv[data.time_Gyr>1.5]
    mask_conv = ~ma.array(data.zbot[data.time_Gyr>1.5], mask = data.conv[data.time_Gyr>1.5]).mask
    
    zbot_cond = ma.array(data.zbot[data.time_Gyr>1.5], mask = mask_cond)
    zbot_conv = ma.array(data.zbot[data.time_Gyr>1.5], mask = mask_conv)
    min_zbot.append(np.min(zbot_conv))
    max_zbot.append(np.max(zbot_conv))

model_list = ['Europa-tidal1_eta1.0d14_P0.1TW_1.5wt%-NH3', # power
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3',
              'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta1.0d14_P1.2TW_1.5wt%-NH3',]
model_list = ['Europa-tidal1_eta3.2d14_P0.6TW_1.5wt%-NH3', # power
              'Europa-tidal1_eta3.2d14_P0.8TW_1.5wt%-NH3',
              'Europa-tidal1_eta3.2d14_P1.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta3.2d14_P1.2TW_1.5wt%-NH3',]


for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax2.plot(x,zbot_conv,color=colors[i],lw=3)
    ax2.plot(x,zbot_cond,color=colors[i],linestyle='dashed')  
#  ------------------------------ figure setting ------------------------------
ax.set_ylim(161,0)
ax2.set_ylim(161,0)
ax2.legend(fontsize=labelsize)
ax2.set_xlabel('time (Gyr)',fontsize=labelsize)

for aa in [ax,ax2]:
    aa.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
playerA = np.array(max_zbot)
playerB = np.array(min_zbot)
hat_graph(axx3, label_list, [playerA, playerB], ['Player A', 'Player B'])
axx3.set_ylim(161,0)
axx3.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
axx3.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=10)
for aa in [axx,axx3]:
    axx3.set_xlabel('internal power (TW)',fontsize=labelsize)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    aa.grid()
fig.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figureS3_v2.pdf')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 20:04:12 2023

@author: chingchen
"""


import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt

labelsize = 30
bwith = 3
variation_eta = 0
variation_P = 1
variation_P_fixvol = 0
pmelt_plot = 0
#----------------------------------------------------------------------------------------
pcrit = 0
pmelt = 0
variation_thermal_evolution = 0
Hvar_tprofile = 0
Hvar_shell_propertise = 0
Hvar_thermal_evolution = 1
plot_final=0

### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
    modelpath = '/Users/chingchen/Desktop/model/'
    figpath = '/Users/chingchen/Desktop/figure/'
else:
    path = '/lfs/jiching/data/'
    workpath = '/lfs/jiching/ScalingLaw_model/'
    modelpath = '/lfs/jiching/ScalingLaw_model/'
    figpath = '/lfs/jiching/figure/'

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']

colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']

if variation_eta:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list = ['Europa-tidal1_eta1.0d12-1.0d15_P0.1TW_1.0wt%-NH3',
                  'Europa-tidal1_eta1.0d12-1.0d15_P0.3TW_1.0wt%-NH3',
                  'Europa-tidal1_eta1.0d12-1.0d15_P0.5TW_1.0wt%-NH3',
                  'Europa-tidal1_eta1.0d12-1.0d15_P0.6TW_1.0wt%-NH3',
                  'Europa-tidal1_eta1.0d12-1.0d15_P0.8TW_1.0wt%-NH3',
                  'Europa-tidal1_eta1.0d12-1.0d15_P1.0TW_1.0wt%-NH3',]
                  # 'Europa-tidal1_eta1d12-1d15_P1.2TW_2.5wt%-NH3',
                  # 'Europa-tidal1_eta1d12-1d15_P1.4TW_2.5wt%-NH3']
    # model_list = ['Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.1TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.3TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.5TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.8TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P1.0TW_1.0wt%-NH3',]
    # model_list = ['Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_0.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_0.5wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_1.5wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_2.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_2.5wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta1d12-1d15_P0.6TW_3.0wt%-NH3',]
    rainbow = cm.get_cmap('winter',len(model_list))
    newcolors = rainbow(np.linspace(0, 1, len(model_list)))
    label_list=['P = 0.1 Tw','P = 0.3 Tw','P = 0.5 Tw','P = 0.6 Tw','P = 0.8 Tw','P = 1.0 Tw']
    # label_list=['vol = 0.0%','vol = 0.5%','vol = 1.0%','vol = 1.5%','vol = 2.0%','vol = 2.5%','vol = 3.0%']
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        eta=data.etaref
        x = eta
        ax.plot(x[data.conv==1],data.zbot[data.conv==1],color=newcolors[i],label = label_list[i],lw=4)
        ax2.plot(x[data.conv==1],data.Tm[data.conv==1],color=newcolors[i],lw=4)
        ax3.plot(x[data.conv==1],data.dlid[data.conv==1],color=newcolors[i],lw=4)
        # ax.plot(x[data.conv==1],data.Ftop[data.conv==1],color=newcolors[i],label = label_list[i],lw=4)
        # ax2.plot(x[data.conv==1],data.Tbot[data.conv==1],color=newcolors[i],lw=4)
        # ax3.plot(x[data.conv==1],data.T_core[data.conv==1],color=newcolors[i],lw=4)
        
        if len(data.zbot[data.conv==0])>0:
            ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            # ax.plot(x[data.conv==0],data.Ftop[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            # ax2.plot(x[data.conv==0],data.Tbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            # ax3.plot(x[data.conv==0],data.T_core[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
            
            ax2.vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                        ymax=data.Tm[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
            ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                        ymax=data.dlid[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    # ax.set_ylabel('bottom heat flux',fontsize = labelsize-5)
    # ax2.set_ylabel('bottom temperature (K)',fontsize = labelsize-5)
    # ax3.set_ylabel(' temperature of core (K)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    ax2.set_ylim(180,280)
    ax3.set_ylim(0,35)
    ax.legend(fontsize = 20)
    ax3.set_xlabel('log(reference viscosity)',fontsize=labelsize)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(1e12,1e15)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if variation_P:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','vol',
                  'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list = ['Europa-tidal1_eta1.0d12_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta3.2d12_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta1.0d13_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta3.2d13_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta1.0d14_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta3.2d14_P0.0-3.0TW_2.5wt%-NH3',
                 'Europa-tidal1_eta1.0d15_P0.0-3.0TW_2.5wt%-NH3',]
    # model_list = ['Europa-tidal1_eta1.0d12_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta3.2d12_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta1.0d13_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta3.2d13_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta1.0d14_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta3.2d14_P0.0-1.5TW_3.0wt%-NH3',
    #              'Europa-tidal1_eta1.0d15_P0.0-1.5TW_3.0wt%-NH3',]
    # model_list = ['Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_0.5wt%-NH3',
    #              'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_1.0wt%-NH3',
    #              'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_1.5wt%-NH3',
    #              'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_2.0wt%-NH3',
    #              'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_2.5wt%-NH3',
    #              'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_3.0wt%-NH3']
   # model_list = ['Europa-tidal5_period0.15Gyr_eta1.0d12_P0-1.5TW_1.5wt%-NH3',
   #               'Europa-tidal5_period0.15Gyr_eta3.2d12_P0-1.5TW_1.5wt%-NH3',
   #               'Europa-tidal5_period0.15Gyr_eta1.0d13_P0-1.5TW_1.5wt%-NH3',
   #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_1.5wt%-NH3',
   #               'Europa-tidal5_period0.15Gyr_eta1.0d14_P0-1.5TW_1.5wt%-NH3',
   #               'Europa-tidal5_period0.15Gyr_eta3.2d14_P0-1.5TW_1.5wt%-NH3']
    rainbow = cm.get_cmap('winter',len(model_list))
    newcolors = rainbow(np.linspace(0, 1, len(model_list)))
   
    label_list=['eta = 1e12','eta = 3.2e12','eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14','eta = 1e15']
    # label_list=['vol = 0.0%','vol = 0.5%','vol = 1.0%','vol = 1.5%','vol = 2.0%','vol = 2.5%','vol = 3.0%']
    for i, model in enumerate(model_list):
       data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
           header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
       data = data.replace('D','e',regex=True).astype(float) # data convert to float
       Pint=data.Pint/1e12
       x = Pint
       ax.plot(x[data.conv==1],data.zbot[data.conv==1],c=newcolors[i],label = label_list[i],lw=5)
       ax2.plot(x[data.conv==1],data.Tm[data.conv==1],color=newcolors[i],lw=5)
       ax3.plot(x[data.conv==1],data.dlid[data.conv==1],color=newcolors[i],lw=5)
       
       ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
       ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
       ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
       
       if len(x[data.conv==0])!=0:
           xx=np.array(x[data.conv==0])[0]
           
           ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                      ymax=data.zbot[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
           ax2.vlines(x=xx, ymin=data.Tm[data.conv==0].iloc[0], 
                      ymax=data.Tm[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
           ax3.vlines(x=xx, ymin=data.dlid[data.conv==0].iloc[0], 
                      ymax=data.dlid[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    ax2.set_ylim(180,300)
    ax3.set_ylim(0,30)
    # ax.legend(fontsize = 20)
    ax3.set_xlabel('P$_{tide}$',fontsize=labelsize)
    ax4  = fig.add_axes([0.95,0.12,0.03,0.76])
    norm = mpl.colors.Normalize(vmin=0,vmax=3)
    cb1  = mpl.colorbar.ColorbarBase(ax4,cmap=rainbow,norm=norm,orientation='vertical')
    cb1.set_label('vol %',fontsize=labelsize)
    cb1.ax.tick_params(axis='y', labelsize=labelsize-2)
    cb1.ax.yaxis.set_label_position('right')
    for aa in [ax,ax2,ax3]:
       aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
       aa.set_xlim(0.0,1.5)
       aa.grid()
       for axis in ['top','bottom','left','right']:
           aa.spines[axis].set_linewidth(bwith)

if variation_P_fixvol:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','vol',
                  'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
   
    model_list = ['Europa-tidal5_period0.15Gyr_eta1.0d12_P0-1.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d12_P0-1.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta1.0d13_P0-1.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0-1.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta1.0d14_P0-1.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d14_P0-1.5TW_1.5wt%-NH3']
    model_list = ['Europa-tidal1_eta1.0d12_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta3.2d12_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta1.0d13_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta3.2d13_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta1.0d14_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta3.2d14_P0.0-1.5TW_3.0wt%-NH3',
                 'Europa-tidal1_eta1.0d15_P0.0-1.5TW_3.0wt%-NH3',]
    
    label_list=['eta = 1e12','eta = 3.2e12','eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14','eta = 1e15']
    for i, model in enumerate(model_list):
       data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
           header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
       data = data.replace('D','e',regex=True).astype(float) # data convert to float
       Pint=data.Pint/1e12
       x = Pint
       qqq=ax.plot(x[data.conv==1],data.zbot[data.conv==1],c=colors[i],label = label_list[i],lw=5)
       ax2.plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=5)
       ax3.plot(x[data.conv==1],data.dlid[data.conv==1],color=colors[i],lw=5)
       
       ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       
       if len(x[data.conv==0])!=0:
           xx=np.array(x[data.conv==0])[0]
           
           ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                      ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
           ax2.vlines(x=xx, ymin=data.Tm[data.conv==0].iloc[0], 
                      ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
           ax3.vlines(x=xx, ymin=data.dlid[data.conv==0].iloc[0], 
                      ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    ax2.set_ylim(180,300)
    ax3.set_ylim(0,30)
    ax.legend(fontsize = 15)
    ax3.set_xlabel('P$_{tide}$',fontsize=labelsize)
    for aa in [ax,ax2,ax3]:
       aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
       aa.set_xlim(0.0,1.5)
       aa.grid()
       for axis in ['top','bottom','left','right']:
           aa.spines[axis].set_linewidth(bwith)

if pmelt_plot: #第二組不同厚度ice的結果皆相同 不能用這種畫法
    header_list = ['Pint','Hint','etaref','pCLA','conv','zbot','%vol',
                   'Tbot','Tm','Fbot','Ftop','dlid']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model = 'Europa-tidal1_eta1e13_P0.0-1.5TW_3.0wt%-NH3'
    
    data = pd.read_csv(workpath+model+'_Pcrit-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12

    ax.scatter(Pint, data.zbot,color=newcolors[5],label = 'eta = 1e13')
    ax2.scatter(Pint,data.Tm,color=newcolors[5])
    ax3.scatter(Pint, data.dlid,color=newcolors[5])
    
    model = 'Europa_eta1e14_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_Pcrit-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.zbot,color=newcolors[10],label = 'eta = 1e14')
    ax2.scatter(Pint,data.Tm,color=newcolors[10])
    ax3.scatter(Pint, data.dlid,color=newcolors[10])
    
    model = 'Europa_eta1e15_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_Pcrit-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.zbot,color=newcolors[2],label = 'eta = 1e15')
    ax2.scatter(Pint,data.Tm,color=newcolors[2])
    ax3.scatter(Pint, data.dlid,color=newcolors[2])
    
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax3.set_xlabel('P$_{tide}$',fontsize=labelsize)
        #aa.set_xlim(0.1,3)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
                
    fig2,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model = 'Europa_eta1e14_P0.01-3TW_D40_3.0%-NH3_k2.6_wt%'
    
    data = pd.read_csv(workpath+model+'_Pmelt-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.Ftop,color=newcolors[5],label = 'D = 40 km')
    ax2.scatter(Pint,data.Tm,color=newcolors[5])
    ax3.scatter(Pint, data.dlid,color=newcolors[5])
    
    model = 'Europa_eta1e14_P0.01-3TW_D80_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_Pmelt-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.Ftop,color=newcolors[10],label = 'D = 80 km')
    ax2.scatter(Pint,data.Tm,color=newcolors[10])
    ax3.scatter(Pint, data.dlid,color=newcolors[10])
    
    model = 'Europa_eta1e14_P0.01-3TW_D120_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_Pmelt-vs-zbot.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.Ftop,color=newcolors[2],label = 'D = 120 km')
    ax2.scatter(Pint,data.Tm,color=newcolors[2])
    ax3.scatter(Pint, data.dlid,color=newcolors[2])
    
    ax.set_ylabel('surface heat flux',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    #ax.set_ylim(125,0)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax3.set_xlabel('Pint',fontsize=labelsize)
        aa.grid()
        #aa.set_xlim(0.003,3)
        #aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if pcrit:
    header_list = ['Pint','Hint','etaref','pCLA','conv','zbot','%vol',
                   'Tbot','Tm','Fbot','Ftop','dlid']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list=['Europa-tidal1_eta1e13_P0.0-1.5TW_3.0wt%-NH3',
                ]
    
    label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0 Tw','P = 1 Tw','P = 2 Tw','P = 3 Tw',]
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Pcrit-vs-zbot.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint
        x = Pint
        # x=data.pCLA * 100
        ax.scatter(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i])
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i])
        ax3.scatter(x[data.conv==1], data.dlid[data.conv==1],color=colors[i])
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    # ax.set_ylim(161,0)
    # ax2.set_ylim(230,300)
    # ax3.set_ylim(0,30)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax3.set_xlabel('pCLA',fontsize=labelsize)
        # aa.set_xlim(0.01,3)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        aa.grid()
        # aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if pmelt:
    header_list = ['Pint','Hint','etaref','pCLA','conv','zbot','%vol',
                   'Tbot','Tm','Fbot','Ftop','dlid']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list=['Europa-tidal1_eta1e13_P0.0-1.5TW_3.0wt%-NH3',
                ]
    
    label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0 Tw','P = 1 Tw','P = 2 Tw','P = 3 Tw',]
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Pcrit-vs-zbot.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint
        x = Pint
        # x=data.pCLA * 100
        ax.scatter(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i])
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i])
        ax3.scatter(x[data.conv==1], data.dlid[data.conv==1],color=colors[i])
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    # ax.set_ylim(161,0)
    # ax2.set_ylim(230,300)
    # ax3.set_ylim(0,30)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax3.set_xlabel('pCLA',fontsize=labelsize)
        # aa.set_xlim(0.01,3)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        aa.grid()
        # aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)

if variation_thermal_evolution: # 同 variation_P or variation_vis
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list=['Europa-tidal1_eta1e13_P0.0-1.5TW_3.0wt%-NH3',
                ]
    
    label_list=['eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    # label_list=['P = 0.3 Tw','P = 0.6 Tw','P = 1.0 Tw','P = 3 Tw',]
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint/1e12
        x = Pint
        # x=data.pCLA * 100
        # x=data.etaref
        ax.scatter(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i])
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i])
        ax3.scatter(x[data.conv==1], data.dlid[data.conv==1],color=colors[i])
        ax.plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3.plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        
        # ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
        #            ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
        # ax2.vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
        #            ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
        # ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
        #            ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    # ax2.set_ylim(230,300)
    ax3.set_ylim(0,30)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax3.set_xlabel('eta',fontsize=labelsize)
        # aa.set_xlim(0.01,3)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        aa.grid()
        # aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if Hvar_tprofile:
    model = 'Europa-tidal1_eta1e13_P0.3TW_3.0wt%-NH3'
    
    line = open(workpath+model+'_Hvar_T-profiles.dat').readlines(0)
    time = line[1].split("at ")[1].split('\n')[0].split('   ')[1:]
    data = np.loadtxt(workpath+model+'_Hvar_T-profiles.dat',skiprows=2)
    kk = np.hsplit(data, 31)
    
    for ii, mm in enumerate(range(len(kk))):
        uu=pd.DataFrame(kk[ii], columns = ['radius','temperature'])
        uu.to_csv(path+model+'_Tprofiles_'+time[ii], sep=',', index=False)
        
        dd = pd.read_csv(path+model+'_Tprofiles_'+time[ii])
        fig,(ax) = plt.subplots(1,1,figsize=(8,12))
        
        ax.scatter(dd.temperature, dd.radius)
        ax.set_ylim(1300,1561)
        for aa in [ax]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True,top=True,direction='in',pad=15)
            #aa.set_xlim(xmin,xmax)
            #aa.set_ylim(-zmin,-zmax)
            #aa.set_xlim(0,5)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)

if Hvar_shell_propertise: 
    header_list = ['zbot','P_MPa','Hint','conv','melt','pCLA','Tbot','Tm',
                   'etam','k','Fcond_bot','Fcond_top','Fbot','Ftop','dlid',
                   'z_topTBL','z_botTBL']
    model_list = ['Europa-tidal1_eta1e13_P0.3TW_3.0wt%-NH3']
    fig,(ax) = plt.subplots(1,1,figsize=(8,12))
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_shell-properties.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        ax.scatter(data.Tm, data.zbot)
        for aa in [ax]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True,top=True,direction='in',pad=15)
            #aa.set_xlim(xmin,xmax)
            aa.set_ylim(0,161)
            #aa.set_xlim(0,5)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)


if Hvar_thermal_evolution:
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

    # model_list = ['Europa-tidal5_period0.15Gyr_eta3.2d13_P0.1TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.3TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.5TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.6TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.8TW_1.0wt%-NH3',
    #               'Europa-tidal5_period0.15Gyr_eta3.2d13_P1.0TW_1.0wt%-NH3',]
    
    
    model_list = ['Europa-tidal5_period0.15Gyr_eta3.2d13_P0.1TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.3TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P1.0TW_1.5wt%-NH3',]
    
    model_list = ['Europa-tidal5_period0.15Gyr_eta1.0d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta1.0d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta1.0d14_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.15Gyr_eta3.2d14_P0.6TW_1.5wt%-NH3',]
    # Europa-tidal5_period0.15Gyr_etaXXX_P0.6TW_1.5wt%-NH3
    
    model_list = ['Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d14_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d14_P0.6TW_1.5wt%-NH3',]
    
    model_list = ['Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d12_P1.0TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d12_P1.0TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d13_P1.0TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d13_P1.0TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d14_P1.0TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.2_eta3.2d14_P1.0TW_1.5wt%-NH3',]
    
    
    model_list = ['Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d12_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d13_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d14_P0.6TW_1.5wt%-NH3',]
    
    model_list = ['Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d12_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d12_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d13_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d13_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta3.2d14_P0.8TW_1.5wt%-NH3',]
    
    # model_list = ['Europa-tidal5_period0.15Gyr_eta1.0d14_P0.6TW_1.5wt%-NH3',
    #               'Europa-tidal5_period0.14Gyr_emx1.2_eta1.0d14_P0.6TW_1.5wt%-NH3',
    #               'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.6TW_1.5wt%-NH3',]
    
    model_list = ['Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.1TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.3TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.5TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.6TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P0.8TW_1.5wt%-NH3',
                  'Europa-tidal5_period0.14Gyr_emx1.5_eta1.0d14_P1.0TW_1.5wt%-NH3',]
    
    model_list = ['Europa-tidal1_eta3.2d12_P0.3TW_2.5wt%-NH3',
                  'Europa-tidal1_eta3.2d12_P0.6TW_2.5wt%-NH3']
    
    label_list=['eta = 1e12','eta = 3.2e12','eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    # label_list=['eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0.1 Tw','P = 0.3 Tw','P = 0.5 Tw','P = 0.6 Tw','P = 0.8 Tw','P = 1.0 Tw']
    # label_list = ['emx = 1.1','emx = 1.2','emx = 1.5','1.0 Gyr']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(20,18))
    for i, model in enumerate(model_list):
        i = i
        # data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
        #     header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        
        zbot_cond = ma.array(data.zbot, mask = mask_cond)
        zbot_conv = ma.array(data.zbot, mask = mask_conv)
        tm_cond   = ma.array(data.Tm, mask = mask_cond)
        tm_conv   = ma.array(data.Tm, mask = mask_conv)
        dlid_cond = ma.array(data.dlid, mask = mask_cond)
        dlid_conv = ma.array(data.dlid, mask = mask_conv)
        
        
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
        ax2.plot(x,tm_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.plot(x,dlid_conv,color=colors[i],label=label_list[i],lw=3)
        # ax3.axhline(y=6.3, color='gray', linestyle='--',lw = 2)
        
        ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
        ax2.plot(x,tm_cond,color=colors[i],linestyle='dashed')
        ax3.plot(x,dlid_cond,color=colors[i],linestyle='dashed')
        
     
        # ax.plot(x,data.Pint,color=colors[i],lw=3)
        
        ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
        ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
        ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
        ax.set_ylim(161,0)
        ax2.set_ylim(200,300)
        ax3.set_ylim(30,0)
        ax3.legend(fontsize=labelsize)
        for aa in [ax,ax2,ax3]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
            ax3.set_xlabel('Time',fontsize=labelsize)
            aa.set_xlim(0,4.55)
            aa.grid()
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)



if plot_final:
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(24,12))
    for kk,vol in enumerate(['0.5','1.0','1.5','2.0','2.5','3.0',]):
        for power in (['0.1','0.3','0.5','0.6','0.8','1.0',]):
            for vis in (['1.0d12','3.2d12','1.0d13','3.2d13','1.0d14','3.2d14']):
                model  = 'Europa-tidal1_eta'+vis+'_P'+power+'TW_'+vol+'wt%-NH3'
                data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
                    header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True).tail(1) 
                data = data.replace('D','e',regex=True).astype(float) # data convert to float
                mm = int(float(vis.replace('d','e')))
                if len(data.T_core)==0:
                    print(model)
                    continue
                ax.scatter(mm,data.T_core,color=colors[kk],s=100)
                ax2.scatter(mm,data.Tbot,color=colors[kk],s=100)
                ax3.scatter(mm,data.zbot,color=colors[kk],s=100)
        ax.set_xlim(0.5e12,1e15)
    
    ax.set_title('final core temperature',fontsize=labelsize)
    ax.set_xlabel('viscosity',fontsize=labelsize)
    ax.set_ylabel('temperature (K)',fontsize=labelsize-10)
    ax2.set_title('final temperature at bottom ice shell',fontsize=labelsize-8)
    ax2.set_xlabel('viscosity',fontsize=labelsize)
    ax2.set_ylabel('temperature (K) at bottom',fontsize=labelsize-10)
    ax3.set_title('final depth of ice shell',fontsize=labelsize)
    ax3.set_xlabel('viscosity',fontsize=labelsize)
    ax3.set_ylabel('depth (km)',fontsize=labelsize-10)
    for aa in [ax,ax2,ax3]:
        aa.grid()
        aa.set_xscale('log')
        aa.tick_params(labelsize=labelsize-10,width=3,length=10,right=True, top=True,direction='in',pad=15)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
            
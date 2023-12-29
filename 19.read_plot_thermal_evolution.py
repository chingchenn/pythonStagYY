#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 20:04:12 2023

@author: chingchen
"""


import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt

labelsize = 30
bwith = 3
variation_eta = 0
variation_P = 0
pmelt_plot = 0
#----------------------------------------------------------------------------------------
pcrit = 0
pmelt = 0
variation_thermal_evolution = 0
tprofile = 1
shell_propertise = 0
variation_time_thermal_evolution = 0

### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
    # workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution_vs_H-and-NH3/'
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

colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67']

if variation_eta:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list = ['Europa_eta1e12-1e15_P0.0TW_D40_3.0%-NH3_k2.6_wt%',
                  #'Europa_eta1e12-1e15_P0.2TW_D40_3.0%-NH3_k2.6_wt%',
                  'Europa_eta1e12-1e15_P0.3TW_D40_3.0%-NH3_k2.6_wt%',
                  'Europa_eta1e12-1e15_P0.6TW_D40_3.0%-NH3_k2.6_wt%',
                  'Europa_eta1e12-1e15_P1.0TW_D40_3.0%-NH3_k2.6_wt%',
                  'Europa_eta1e12-1e15_P1.2TW_D40_3.0%-NH3_k2.6_wt%']
   
    label_list=['P = 0 Tw','P = 0.2 Tw','P = 0.3 Tw','P = 0.6 Tw','P = 1.0 Tw','P = 1.2 Tw']
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        eta=data.etaref
        x = eta
        ax.plot(x[data.conv==1],data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=4)
        ax2.plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=4)
        ax3.plot(x[data.conv==1],data.dlid[data.conv==1],color=colors[i],lw=4)
        
        ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        if len(data.zbot[data.conv==0])>0:
            ax2.vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                        ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                        ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax.set_ylim(161,0)
    ax2.set_ylim(180,280)
    ax3.set_ylim(0,35)
    ax.legend(fontsize = 28)
    ax3.set_xlabel('log(reference viscosity)',fontsize=labelsize)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(1e12,1e15)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if variation_P:
   header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                  'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
   fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
   model_list = ['Europa_eta3e13_P0-2.0TW_D40_3.0%-NH3_k2.6_wt%',
                 'Europa_eta1e14_P0-2.0TW_D40_3.0%-NH3_k2.6_wt%',
                 'Europa_eta3e14_P0-2.0TW_D40_3.0%-NH3_k2.6_wt%']
  
   
   label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
   for i, model in enumerate(model_list):
       data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
           header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
       data = data.replace('D','e',regex=True).astype(float) # data convert to float
       Pint=data.Pint/1e12
       x = Pint
       ax.plot(x[data.conv==1],data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=5)
       ax2.plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=5)
       ax3.plot(x[data.conv==1],data.dlid[data.conv==1],color=colors[i],lw=5)
       
       ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
       
       
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
   ax.legend(fontsize = 28)
   ax3.set_xlabel('P$_{tide}$',fontsize=labelsize)
   for aa in [ax,ax2,ax3]:
       aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
       aa.set_xlim(0.0,2)
       aa.grid()
       
       for axis in ['top','bottom','left','right']:
           aa.spines[axis].set_linewidth(bwith)


if pmelt_plot: #第二組不同厚度ice的結果皆相同 不能用這種畫法
    header_list = ['Pint','Hint','etaref','pCLA','conv','zbot','%vol',
                   'Tbot','Tm','Fbot','Ftop','dlid']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model = 'Europa_eta1e13_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    
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
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
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
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
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

if variation_thermal_evolution:    
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
                ]
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.3TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.0TW_3.0%-NH3_k2.6_wt%']
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e13_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta3e13_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e14_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',]
    
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
        
        
        ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                   ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
        ax2.vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                   ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
        ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                   ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
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
if tprofile:
    model = 'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%'
    
    
    line = open(workpath+model+'_Hvar_T-profiles.dat').readlines(0)
    #line = file.readlines()
    time = line[1].split("at ")[1].split('\n')[0].split('   ')[1:]
    data = np.loadtxt(workpath+model+'_Hvar_T-profiles.dat',skiprows=2)
    kk = np.hsplit(data, 22)
    
    for ii, mm in enumerate(range(len(kk))):
        uu=pd.DataFrame(kk[ii], columns = ['radius','temperature'])
        uu.to_csv(path+model+'_Tprofiles_'+time[ii], sep=',', index=False)
        
        dd = pd.read_csv(path+model+'_Tprofiles_'+time[ii])
        fig,(ax) = plt.subplots(1,1,figsize=(8,12))
        
        ax.scatter(dd.temperature, dd.radius)
        ax.set_ylim(1300,1561)
        for aa in [ax]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
            #aa.set_xlim(xmin,xmax)
            #aa.set_ylim(-zmin,-zmax)
            #aa.set_xlim(0,5)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
            
            
            
    # data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
    #     header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    # data = data.replace('D','e',regex=True).astype(float) # data convert to float
    # x = data.time_Gyr
    # # ax.scatter(x[data.conv==1], data.zbot[data.conv==1])
    # # ax.scatter(x[data.conv==1],data.zbot[data.conv==1],color=colors[i])
    # ax2.scatter(x[data.conv==1],data.Tm[data.conv==1])#,label=label_list[i])
    # ax3.scatter(x[data.conv==1],data.dlid[data.conv==1])
    # # ax.scatter(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed')
    # ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
    # ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
    
    
    # ax.plot(x,data.Pint,color=colors[i],lw=3)
    
    # # ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    # ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    # ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
    # # ax.set_ylim(161,0)
    # # ax2.set_ylim(190,300)
    # # ax3.set_ylim(0,30)
    # ax2.legend(fontsize=labelsize)
    # for aa in [ax,ax2,ax3]:
    #     aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    #     ax3.set_xlabel('Time',fontsize=labelsize)
    #     aa.set_xlim(0,4.55)
    #     aa.grid()
    #     for axis in ['top','bottom','left','right']:
    #         aa.spines[axis].set_linewidth(bwith)
shell_propertise = 1
if variation_time_thermal_evolution:
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']
    model_list=['Europa-tidal1_period0.0Gyr_m0.0_core1.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.2_core0.8_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.5_core0.5_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.7_core0.3_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%']
    # model_list = ['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.6TW_3.0%-NH3_k2.6_wt%']
    # model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['x2-Europa_eta1e14_3.0%-NH3_k2.6_P0.3TW_wt%',
    #               'xtest-period1.0Gyr-Europa_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'xtest-period0.25Gyr-Europa_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%']
    # # model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P1.0TW_3.0%-NH3_k2.6_wt%',
    # #               'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P1.0TW_3.0%-NH3_k2.6_wt%',
    # #               'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P1.0TW_3.0%-NH3_k2.6_wt%',
    # #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P1.0TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',]

    model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',]
    
    
    label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0 Tw','P = 1 Tw','P = 2 Tw','P = 3 Tw',]
    label_list = ['0.125 Gyr','0.25 Gyr','0.5 Gyr','1.0 Gyr']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(20,18))
    for i, model in enumerate(model_list):
        i = i
        # data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
        #     header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        # ax.scatter(x[data.conv==1], data.zbot[data.conv==1])
        # ax.scatter(x[data.conv==1],data.zbot[data.conv==1],color=colors[i])
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i])#,label=label_list[i])
        ax3.scatter(x[data.conv==1],data.dlid[data.conv==1],color=colors[i])
        # ax.scatter(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed')
        ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        
        ax.plot(x,data.Pint,color=colors[i],lw=3)
        
        # ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
        ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
        ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize-5)
        # ax.set_ylim(161,0)
        # ax2.set_ylim(190,300)
        # ax3.set_ylim(0,30)
        ax2.legend(fontsize=labelsize)
        for aa in [ax,ax2,ax3]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
            ax3.set_xlabel('Time',fontsize=labelsize)
            aa.set_xlim(0,4.55)
            aa.grid()
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
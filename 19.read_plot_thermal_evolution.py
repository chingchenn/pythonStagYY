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


model = 'Europa_eta1e12-1e15_P0.1TW_3.0%-NH3_k2.6_wt%'
labelsize = 30
bwith = 3
variation_eta = 1
variation_P = 1

### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution_vs_H-and-NH3/'
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

if variation_eta:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    # read evolution.dat file
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    
    ax.scatter(data.etaref, data.zbot,color=newcolors[5],label = 'P = 0.1 TW')
    ax2.scatter(data.etaref,data.Tm,color=newcolors[5])
    ax3.scatter(data.etaref, data.dlid,color=newcolors[5])

    model = 'Europa_eta1e12-1e15_P0.6TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(data.etaref, data.zbot,color=newcolors[10],label = 'P = 0.6 TW')
    ax2.scatter(data.etaref,data.Tm,color=newcolors[10])
    ax3.scatter(data.etaref, data.dlid,color=newcolors[10])
    
    model = 'Europa_eta1e12-1e15_P1.8TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(data.etaref, data.zbot,color=newcolors[2],label = 'P = 1.8 TW')
    ax2.scatter(data.etaref,data.Tm,color=newcolors[2])
    ax3.scatter(data.etaref, data.dlid,color=newcolors[2])
    
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax2.set_xlabel('etaref',fontsize=labelsize)
        #aa.set_xlim(xmin,xmax)
        #aa.set_ylim(-zmin,-zmax)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax.set_ylim(161,0)
if variation_P:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model = 'Europa_eta1e13_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float

    ax.scatter(data.Pint, data.zbot,color=newcolors[5],label = 'eta = 1e13')
    ax2.scatter(data.Pint,data.Tm,color=newcolors[5])
    ax3.scatter(data.Pint, data.dlid,color=newcolors[5])
    
    model = 'Europa_eta1e14_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(data.Pint, data.zbot,color=newcolors[10],label = 'eta = 1e14')
    ax2.scatter(data.Pint,data.Tm,color=newcolors[10])
    ax3.scatter(data.Pint, data.dlid,color=newcolors[10])
    
    model = 'Europa_eta1e15_P0.01-3TW_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(data.Pint, data.zbot,color=newcolors[2],label = 'eta = 1e15')
    ax2.scatter(data.Pint,data.Tm,color=newcolors[2])
    ax3.scatter(data.Pint, data.dlid,color=newcolors[2])
    
    ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('interior temperature (K)',fontsize = labelsize)
    ax3.set_ylabel('Stagnant lid thickness (km',fontsize = labelsize)
    ax.set_ylim(161,0)
    ax.legend(fontsize = 18)
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax2.set_xlabel('etaref',fontsize=labelsize)
        aa.set_xlim(1e10,3e12)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
            
            
    fig2,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(12,20))
    model = 'Europa_eta1e14_P0.01-3TW_D40_3.0%-NH3_k2.6_wt%'
    
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    ax.scatter(Pint, data.Ftop,color=newcolors[5],label = 'eta = 1e15')
    ax2.scatter(Pint,data.Tm,color=newcolors[5])
    ax3.scatter(Pint, data.dlid,color=newcolors[5])
    
    model = 'Europa_eta1e14_P0.01-3TW_D80_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(Pint, data.Ftop,color=newcolors[10],label = 'eta = 1e15')
    ax2.scatter(Pint,data.Tm,color=newcolors[10])
    ax3.scatter(Pint, data.dlid,color=newcolors[10])
    
    model = 'Europa_eta1e14_P0.01-3TW_D120_3.0%-NH3_k2.6_wt%'
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    ax.scatter(Pint, data.Ftop,color=newcolors[2],label = 'eta = 1e15')
    ax2.scatter(Pint,data.Tm,color=newcolors[2])
    ax3.scatter(Pint, data.dlid,color=newcolors[2])
    
    for aa in [ax,ax2,ax3]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax2.set_xlabel('etaref',fontsize=labelsize)
        aa.grid()
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
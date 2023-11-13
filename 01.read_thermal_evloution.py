#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:05:06 2023

@author: chingchen
"""


import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt


model = 'test01'
labelsize = 30
bwith = 3
evolution = 0
shell_properties = 1
T_profile = 0
### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/time_dependent_H_thermal_evloution/'
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

if evolution:
    # read evolution.dat file
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    
    fig,(ax,ax2) = plt.subplots(2,1,figsize=(12,12))
    
    ax.scatter(data.dlid, data.Tbot)
    #ax.scatter(data.time_Gyr, data.Ftop/(0.7**2))
    # ax.scatter(data.time_Gyr, data.Power)
    ax2.scatter(data.time_Gyr,data.Tm)
    ax2.scatter(data.time_Gyr,data.T_core)
    
    for aa in [ax,ax2]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax2.set_xlabel('Time',fontsize=labelsize)
        #aa.set_xlim(xmin,xmax)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)


if shell_properties:
    # read shell-properties.dat file
    header_list = ['zbot','P','Hint','conv','melt','pCLA',
                   'Tbot','Tm','etam','k','Fcond_bot','Fcond_top',
                   'Fbot','Ftop','dlid','z_topTBL','z_botTBL']
    data = pd.read_csv(workpath+model+'_shell-properties.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    
    fig,(axmm,ax2) = plt.subplots(2,1,figsize=(12,12))
    
    axmm.scatter(data.etam, data.Tm)
    #ax.scatter(data.zbot, data.Ftop/(0.7**2))
    # ax.scatter(data.time_Gyr, data.Power)
    #ax2.scatter(data.zbot,data.P)
    #ax2.scatter(data.zbot,data.T_core)
    
    for aa in [axmm,ax2]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax2.set_xlabel('zbot',fontsize=labelsize)
        #aa.set_xlim(xmin,xmax)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)


if T_profile:
    # read T-profiles.dat file
    header_list = ['Radius','temperature']
    line = open(workpath+model+'_T-profiles.dat').readlines(0)
    #line = file.readlines()
    time = line[1].split("at ")[1].split('\n')[0].split('   ')[1:]
    data = np.loadtxt(workpath+model+'_T-profiles.dat',skiprows=2)
    kk = np.hsplit(data, 22)
    
    for ii, mm in enumerate(range(len(kk))):
        uu=pd.DataFrame(kk[ii], columns = ['radius','temperature'])
        uu.to_csv(path+model+'_Tprofiles_'+time[ii], sep=',', index=False)
        
    dd = pd.read_csv(path+model+'_Tprofiles_'+time[1])
    fig,(ax) = plt.subplots(1,1,figsize=(8,12))
    
    ax.scatter(dd.temperature, dd.radius)
    
    for aa in [ax]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        #aa.set_xlim(xmin,xmax)
        #aa.set_ylim(-zmin,-zmax)
        #aa.set_xlim(0,5)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:51:46 2023

@author: chingchen
"""

import stagpy
import numpy as np
from stagpy import field
from stagpy import stagyydata
import matplotlib.pyplot as plt
model = 'h01'
path = '/Users/chingchen/Desktop/data/'
path = '/lfs/jiching/thermo_chemical/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/lfs/jiching/figure/'

#plt.rcParams["font.family"] = "Times New Roman"
data = stagyydata.StagyyData(path+model)
plotting_3field = 0
plotting_Tv = 0
plotting_T = 1
plotting_bs = 0
plotting_prim = 0
gif = 1
mp4 = 1
shot = 1000
if plotting_3field:
   
    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
    eta1,eta2,eta3,eta4 = field.get_meshes_fld(data.snaps[shot],'eta')
    rho1,rho2,rho3,rho4 = field.get_meshes_fld(data.snaps[shot],'rho')

    fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(12,4)) 
    ax.set_aspect('equal')
    cmap = plt.cm.get_cmap('RdBu_r')
    colorbar = ax.pcolormesh(kk1,kk2,kk3,cmap = cmap, vmin = 0,vmax = 1)
    ax.axis('off')
    cax = plt.axes([0.165, 0.05, 0.15, 0.05])
    cc1=fig.colorbar(colorbar, ax=ax,cax=cax,orientation='horizontal')
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='Temperature', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.set_title(model+' at time '+str(shot/1000),fontsize = 26)
    ax2.set_aspect('equal')
    cmap = plt.cm.get_cmap('rainbow')
    colorbar = ax2.pcolormesh(kk1,kk2,np.log10(eta3),cmap = cmap, vmin = -4,vmax = 4)
    ax2.axis('off')
    cax = plt.axes([0.445, 0.05, 0.15, 0.05])
    cc2=fig.colorbar(colorbar, ax=ax2,cax=cax,orientation='horizontal')
    cc2.ax.tick_params(labelsize=20)
    cc2.set_label(label='Viscosity', size=25)
    ax3.set_aspect('equal')
    cmap = plt.cm.get_cmap('afmhot_r')
    colorbar = ax3.pcolormesh(kk1,kk2,rho3,cmap = cmap, vmin = 0.9,vmax = 1.1)
    ax3.axis('off')
    cax = plt.axes([0.72, 0.05, 0.15, 0.05])
    cc3=fig.colorbar(colorbar, ax=ax3,cax=cax,orientation='horizontal')
    cc3.set_label(label='Density', size=25)
    cc3.ax.tick_params(labelsize=20)
    fig.savefig(figpath+model+'_'+'snapshot_'+str(shot)+'_field.png')
    fig.gca()
    plt.close(fig)
if plotting_Tv:
   
    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
    eta1,eta2,eta3,eta4 = field.get_meshes_fld(data.snaps[shot],'eta')

    fig,(ax,ax2) = plt.subplots(1,2,figsize=(8,5)) 
    ax.set_aspect('equal')
    cmap = plt.cm.get_cmap('RdBu_r')
    colorbar = ax.pcolormesh(kk1,kk2,kk3,cmap = cmap, vmin = 0,vmax = 1)
    ax.axis('off')
    cax = plt.axes([0.157, 0.15, 0.3, 0.05])
    cc1=fig.colorbar(colorbar, ax=ax,cax=cax,orientation='horizontal')
    cc1.ax.tick_params(labelsize=15)
    cc1.set_label(label='Temperature', size=15)
    cc1.ax.yaxis.set_label_position('left')
    ax.set_title(model+' at time '+str(shot/1000),fontsize = 26)
    ax2.set_aspect('equal')
    cmap = plt.cm.get_cmap('rainbow')
    colorbar = ax2.pcolormesh(kk1,kk2,np.log10(eta3),cmap = cmap, vmin = -2.5,vmax = 2.5)
    ax2.axis('off')
    cax = plt.axes([0.585, 0.15, 0.30, 0.05])
    cc2=fig.colorbar(colorbar, ax=ax2,cax=cax,orientation='horizontal')
    cc2.ax.tick_params(labelsize=15)
    cc2.set_label(label='Viscosity', size=15)

    fig.savefig(figpath+model+'_'+'snapshot_'+str(shot)+'_field.png')
    fig.gca()
    plt.close(fig)

if plotting_T:
    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
    ### Cause the mesh lack one column, 
    ### need to concatenate it for xmesh, ymseh and field
    kk1 = np.concatenate((kk1, kk1[:1]), axis=0)
    kk2 = np.concatenate((kk2, kk2[:1]), axis=0)
    newline = (kk3[:1] + kk3[-1:]) / 2
    kk3 = np.concatenate((kk3, newline), axis=0)

    # normalized the field to [0-1]
    nor_kk3=(kk3-np.min(kk3))/(np.max(kk3)-np.min(kk3))

    ### plot begin
    fig,(ax) = plt.subplots(1,figsize=(13,13)) 
    ax.set_aspect('equal')
    cmap = plt.cm.get_cmap('RdBu_r')
    colorbar = ax.pcolormesh(kk1,kk2,nor_kk3,cmap = cmap, vmin = 0,vmax = 1)
    ax.axis('off')
    ax.set_title(model+' at time '+str(shot/1000),fontsize = 12)
    cax = plt.axes([0.93, 0.285, 0.01, 0.431])
    cbar = plt.colorbar(colorbar, cax=cax)
    cbar.set_label(label = 'Temperature',size=12)

    fig.savefig(figpath+model+'_'+'Temperature_snapshot_'+str(shot)+'_field.png')
    print('save figure'+model+'_'+'Temperature_snapshot_'+str(shot)+'_field.png')
    fig.gca()
    plt.close(fig)

if plotting_bs:
    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'bs')
    ### Cause the mesh lack one column, 
    ### need to concatenate it for xmesh, ymseh and field
    kk1 = np.concatenate((kk1, kk1[:1]), axis=0)
    kk2 = np.concatenate((kk2, kk2[:1]), axis=0)
    newline = (kk3[:1] + kk3[-1:]) / 2
    kk3 = np.concatenate((kk3, newline), axis=0)

    # normalized the field to [0-1]
    nor_kk3=(kk3-np.min(kk3))/(np.max(kk3)-np.min(kk3))

    ### plot begin
    fig,(ax) = plt.subplots(1,figsize=(13,13)) 
    ax.set_aspect('equal')
    cmap = plt.cm.get_cmap('RdBu_r')
    colorbar = ax.pcolormesh(kk1,kk2,nor_kk3,cmap = cmap, vmin = 0,vmax = 1)
    ax.axis('off')
    ax.set_title(model+' at time '+str(shot/1000),fontsize = 12)
    cax = plt.axes([0.93, 0.285, 0.01, 0.431])
    cbar = plt.colorbar(colorbar, cax=cax)
    cbar.set_label(label = 'Basalt',size=12)

    fig.savefig(figpath+model+'_'+'Basalt_snapshot_'+str(shot)+'_field.png')
    print('save figure'+model+'_'+'Basalt_snapshot_'+str(shot)+'_field.png')
    fig.gca()
    plt.close(fig)

if plotting_prim:
    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'prim')
    ### Cause the mesh lack one column, 
    ### need to concatenate it for xmesh, ymseh and field
    kk1 = np.concatenate((kk1, kk1[:1]), axis=0)
    kk2 = np.concatenate((kk2, kk2[:1]), axis=0)
    newline = (kk3[:1] + kk3[-1:]) / 2
    kk3 = np.concatenate((kk3, newline), axis=0)

    # normalized the field to [0-1]
    nor_kk3=(kk3-np.min(kk3))/(np.max(kk3)-np.min(kk3))

    ### plot begin
    fig,(ax) = plt.subplots(1,figsize=(13,13)) 
    ax.set_aspect('equal')
    cmap = plt.cm.get_cmap('RdBu_r')
    colorbar = ax.pcolormesh(kk1,kk2,nor_kk3,cmap = cmap, vmin = 0,vmax = 1)
    ax.axis('off')
    ax.set_title(model+' at time '+str(shot/1000),fontsize = 12)
    cax = plt.axes([0.93, 0.285, 0.01, 0.431])
    cbar = plt.colorbar(colorbar, cax=cax)
    cbar.set_label(label = 'Primordial Material',size=12)

    fig.savefig(figpath+model+'_'+'Primordial_snapshot_'+str(shot)+'_field.png')
    print('save figure'+model+'_'+'Primordial_snapshot_'+str(shot)+'_field.png')
    fig.gca()
    plt.close(fig)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:54:16 2023

@author: chingchen
"""

import stagpy
import numpy as np
from stagpy import field
from stagpy import stagyydata
import matplotlib.pyplot as plt
model = 'w0801'
path = '/Users/chingchen/Desktop/data/'
path = '/lfs/jiching/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/lfs/jiching/figure/'

plt.rcParams["font.family"] = "Times New Roman"
data = stagyydata.StagyyData(path+model)
plotting_3field = 0
plotting_Tv = 1
gif = 0
mp4 = 0
end = 300
if plotting_3field:
    for shot in range(1,end):
        kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
        eta1,eta2,eta3,eta4 = field.get_meshes_fld(data.snaps[shot],'eta')
        rho1,rho2,rho3,rho4 = field.get_meshes_fld(data.snaps[shot],'rho')

        fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(12,4)) 
        ax.set_aspect('equal')
        cmap = plt.cm.get_cmap('RdBu_r')
        colorbar = ax.scatter(kk1,kk2,c= kk3,cmap = cmap, vmin = 0,vmax = 1)
        ax.axis('off')
        cax = plt.axes([0.165, 0.05, 0.15, 0.05])
        cc1=fig.colorbar(colorbar, ax=ax,cax=cax,orientation='horizontal')
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='Temperature', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax.set_title(model+' at time '+str(shot/1000),fontsize = 26)
        ax2.set_aspect('equal')
        cmap = plt.cm.get_cmap('rainbow')
        colorbar = ax2.scatter(kk1,kk2,c= np.log10(eta3),cmap = cmap, vmin = -4,vmax = 4)
        ax2.axis('off')
        cax = plt.axes([0.445, 0.05, 0.15, 0.05])
        cc2=fig.colorbar(colorbar, ax=ax2,cax=cax,orientation='horizontal')
        cc2.ax.tick_params(labelsize=20)
        cc2.set_label(label='Viscosity', size=25)
        ax3.set_aspect('equal')
        cmap = plt.cm.get_cmap('afmhot_r')
        colorbar = ax3.scatter(kk1,kk2,c= rho3,cmap = cmap, vmin = 0.9,vmax = 1.1)
        ax3.axis('off')
        cax = plt.axes([0.72, 0.05, 0.15, 0.05])
        cc3=fig.colorbar(colorbar, ax=ax3,cax=cax,orientation='horizontal')
        cc3.set_label(label='Density', size=25)
        cc3.ax.tick_params(labelsize=20)
        fig.savefig(figpath+model+'_'+'snapshot_'+str(shot)+'_field.png')
        fig.gca()
        plt.close(fig)
if plotting_Tv:
    print('plotting')
    for shot in range(end-2,end):
        kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
        eta1,eta2,eta3,eta4 = field.get_meshes_fld(data.snaps[shot],'eta')

        fig,(ax,ax2) = plt.subplots(1,2,figsize=(8,5)) 
        ax.set_aspect('equal')
        cmap = plt.cm.get_cmap('RdBu_r')
        colorbar = ax.scatter(kk1,kk2,c= kk3,cmap = cmap, vmin = 0,vmax = 1)
        ax.axis('off')
        cax = plt.axes([0.157, 0.15, 0.3, 0.05])
        cc1=fig.colorbar(colorbar, ax=ax,cax=cax,orientation='horizontal')
        cc1.ax.tick_params(labelsize=15)
        cc1.set_label(label='Temperature', size=15)
        cc1.ax.yaxis.set_label_position('left')
        ax.set_title(model+' at time '+str(shot/1000),fontsize = 26)
        ax2.set_aspect('equal')
        cmap = plt.cm.get_cmap('rainbow')
        colorbar = ax2.scatter(kk1,kk2,c= np.log10(eta3),cmap = cmap, vmin = -2.5,vmax = 2.5)
        ax2.axis('off')
        cax = plt.axes([0.585, 0.15, 0.30, 0.05])
        cc2=fig.colorbar(colorbar, ax=ax2,cax=cax,orientation='horizontal')
        cc2.ax.tick_params(labelsize=15)
        cc2.set_label(label='Viscosity', size=15)

        fig.savefig(figpath+model+'_'+'snapshot_'+str(shot)+'_field.png')
        fig.gca()
        plt.close(fig)
#-----------------------------creat GIF-----------------------------------------
if gif: 
    from PIL import Image
     
    # Create the frames
    frames = []
    for shot in  range(1,end):
        img=figpath+model+'_'+'snapshot_'+str(shot)+'_field.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(figpath+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=40, loop=0)
    
#-----------------------------creat mp4-----------------------------------------    
if mp4:
    import moviepy.editor as mp
    clip = mp.VideoFileClip(figpath+'png_to_gif.gif')
    clip.write_videofile(figpath+'fieldStagYY'+model+".mp4")

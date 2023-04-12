#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 00:10:52 2023

@author: chingchen
"""

import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

model = 'w0207'
path = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/StagYY/'
mp4 = 0
labelsize = 30
bwith = 3
header_list = ['r','Tmean','Tmin','Tmax','vrms','vmin','vmax',
               'vzabs','vzmin','vzmax','vhrms','vhmin','vhmax',
               'etalog','etamin','etamax','elog','emin','emax',
               'slog','smin','smax','whrms','whmin','whmax',
               'wzrms','wzmin','wzmax','drms','dmin','dmax',
               'enadv','endiff','enradh','enviscdiss','enadiabh',
               'cmean','cmin','cmax','rhomean','rhomin','rhomax',
               'airmean','airmin','airmax','primmean','primmin',
               'primmax','ccmean','ccmin','ccmax','fmeltmean',
               'fmeltmin','fmeltmax','metalmean','metalmin',
               'metalmax','gsmean','gsmin','gsmax','viscdisslog',
               'viscdissmin','viscdissmax','advtot','advdesc',
               'advasc','tcondmean','tcondmin','tcondmax']

for i in range(1,2):
    fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(18,15))
    model = 'w0201'
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax.plot(ff.Tmean,ff.r,color = '#414F67',lw=5)
    
    model = 'w0204'
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax2.plot(ff.Tmean,ff.r,color = '#414F67',lw=5)
    
    model = 'w0207'
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax3.plot(ff.Tmean,ff.r,color = '#414F67',lw=5)
    
    for aa in [ax,ax2,ax3]:
        aa.tick_params(axis='x', labelsize=labelsize)
        aa.tick_params(axis='y', labelsize=labelsize)
        aa.spines['bottom'].set_linewidth(bwith)
        aa.spines['top'].set_linewidth(bwith)
        aa.spines['right'].set_linewidth(bwith)
        aa.spines['left'].set_linewidth(bwith)
        aa.grid()
        aa.set_ylim(0,1)
        aa.set_xlim(0,1)
        aa.set_xlabel('Temperature',fontsize = labelsize)
    # ax2.set_xlabel('Temperature',fontsize = labelsize)
    ax.set_ylabel('Depth',fontsize = labelsize)
    
    ax.set_title('Time = '+str(i),fontsize = labelsize)
    ax.text(0.45,0.9,'Ea = 10^5',fontsize = labelsize+8)
    ax2.text(0.45,0.9,'Ea = 10^6',fontsize = labelsize+8)
    ax3.text(0.45,0.9,'Ea = 10^7',fontsize = labelsize+8)
    # fig.savefig(figpath+'frame_'+str(i)+'temperature_profile.png')
    # fig.gca()
    # plt.close(fig)
if mp4 :
    # -----------------------------creat GIF-----------------------------------------
    from PIL import Image
    import glob
     
    # Create the frames
    frames = []
    for i in  range(1,191):
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        
        img=figpath+'frame_'+str(i)+'temperature_profile.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(figpath+'frame_'+qq+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                    save_all=True, duration=60, loop=0)
    
    #-----------------------------creat mp4-----------------------------------------    
    import moviepy.editor as mp
    clip = mp.VideoFileClip(figpath+'frame_'+qq+'png_to_gif.gif')
    clip.write_videofile(figpath+'temperature_profile_'+model+".mp4")
        
    
# fig,(ax) = plt.subplots(1,1,figsize=(6,12))
# newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
#               '#AE6378','#282130','#7E9680','#24788F',
#               '#849DAB','#EA5E51','#35838D','#4198B9',
#               '#414F67','#97795D','#6B0D47','#A80359','#52254F']

# for kk, model in enumerate(['test04','test01','test03']):
#     i = 750
#     ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
#                       sep = '\\s+',header = None,names = header_list)    
#     ax.plot(ff.vzabs*ff.Tmean,ff.r,color = newcolors[kk],lw=5,label = model)
#     # ax.plot(ff.Tmean,ff.r,color = newcolors[kk],lw=5,label = model)
#     ax.tick_params(axis='x', labelsize=labelsize)
#     ax.tick_params(axis='y', labelsize=labelsize)
#     ax.spines['bottom'].set_linewidth(bwith)
#     ax.spines['top'].set_linewidth(bwith)
#     ax.spines['right'].set_linewidth(bwith)
#     ax.spines['left'].set_linewidth(bwith)
#     ax.grid()
#     ax.set_ylim(0,1)
#     # ax.set_xlim(0,1)
#     ax.set_xlabel('Temperature x Vz',fontsize = labelsize)
#     # ax2.set_xlabel('Temperature',fontsize = labelsize)
#     ax.set_ylabel('Depth',fontsize = labelsize)
#     ax.legend(fontsize = 20)
#     ax.set_title('Step = '+str(i),fontsize = labelsize)

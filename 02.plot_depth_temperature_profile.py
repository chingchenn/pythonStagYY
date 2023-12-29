#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 00:10:52 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"


workpath = '/Users/chingchen/Desktop/StagYY_Works/'
path = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/StagYY/temperature_profile/'
datapath = '/Users/chingchen/Desktop/data/'
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
end = 500
#Ra1e4_Ea1e5_f0.75
#Ra1e4_Ea1e6
#Ra1e5_Ea1e6
#Ra1e5_Ea1e7
#Ra3.2e4_Ea1e5
#Ra3.2e4_Ea1e6_f0.7
#Ra3.2e4_Ea1e7_f0.8
#Ra3.2e5_Ea1e5_f0.8


model_information = pd.read_csv(workpath+'model_information_all.csv',sep=',') 
for jj in range(len(model_information)):
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(18,15))
    model = model_information.model[jj]
    end = model_information.time[jj] - 3
    i = end
    rsurf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
#model = 'Ra3.2e4_Ea1e7_f0.8'
    ff = pd.read_csv(path+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax.scatter(ff.Tmean,ff.r,color = '#414F67',s=50)



    ax2.scatter(ff.Tmean,ff.r,color = '#414F67',s=50)
    ax2.axvline(x=Tm,color ='#97795D')
    ax.set_xlim(0,1)
    ax2.set_xlim(0.8,1)
    for aa in [ax,ax2]:
        aa.tick_params( labelsize=labelsize)
        aa.grid()
        aa.set_ylim(0,1)
        
        aa.set_xlabel('Temperature',fontsize = labelsize)
        for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
    #ax2.set_xlabel('Temperature',fontsize = labelsize)
    ax.set_ylabel('Depth',fontsize = labelsize)
    ax.set_title('Rasurf = '+str(rsurf)+' Ea = '+str(format(Ea,'.1e')),fontsize = labelsize)
    ax2.set_title('Tm = '+str(Tm)+'   f = '+str(f),fontsize = labelsize)
    
    # fig.savefig(figpath+model+'frame_'+str(i)+'temperature_profile.png')
    print(figpath+model+'_frame_'+str(i)+'temperature_profile.png')
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
#     ax.tick_params(labelsize=labelsize)
#     for axis in ['top','bottom','left','right']:
#               ax.spines[axis].set_linewidth(bwith)
#     ax.grid()
#     ax.set_ylim(0,1)
#     # ax.set_xlim(0,1)
#     ax.set_xlabel('Temperature x Vz',fontsize = labelsize)
#     # ax2.set_xlabel('Temperature',fontsize = labelsize)
#     ax.set_ylabel('Depth',fontsize = labelsize)
#     ax.legend(fontsize = 20)
#     ax.set_title('Step = '+str(i),fontsize = labelsize)

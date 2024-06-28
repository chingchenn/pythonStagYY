#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 17:44:24 2024

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
plt.rcParams["font.family"] = "Helvetica"
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']

colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']


melting_point_TW=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4])
melting_point_frac=np.array([0.75,0.4,0.3,0.25,0.2,0.16,0.14,0.12,0.11,0.10,0.09,0.08,0.07])
melting_point_depth=[172,243,364,375,381,320,420,370,380,350,380,320,370]


model = 'Europa-tidal1_eta1.0d14_P0.6TW_3.0wt%-NH3'
model = 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3_core0.10_Hvar_2'
model_list = ['Europa-tidal1_eta5.6d13_P1.4TW_1.5wt%-NH3_core0.00_Hvar_2', # power
              'Europa-tidal1_eta5.6d13_P1.4TW_1.5wt%-NH3_core0.02_Hvar_2',
              'Europa-tidal1_eta5.6d13_P1.4TW_1.5wt%-NH3_core0.04_Hvar_2',
              'Europa-tidal1_eta5.6d13_P1.4TW_1.5wt%-NH3_core0.06_Hvar_2',
              'Europa-tidal1_eta5.6d13_P1.4TW_1.5wt%-NH3_core0.07_Hvar_2', # power
              # 'Europa-tidal1_eta5.6d13_P1.2TW_1.5wt%-NH3_core0.09_Hvar_2',
               # 'Europa-tidal1_eta5.6d13_P1.0TW_1.5wt%-NH3_core0.11_Hvar_2',
               # 'Europa-tidal1_eta5.6d13_P0.8TW_1.5wt%-NH3_core0.15_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.7TW_1.5wt%-NH3_core0.16_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.5TW_1.5wt%-NH3_core0.25_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.4TW_1.5wt%-NH3_core0.30_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.3TW_1.5wt%-NH3_core0.35_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.3TW_1.5wt%-NH3_core0.40_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.2TW_1.5wt%-NH3_core0.70_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.2TW_1.5wt%-NH3_core0.75_Hvar_2',
              ]
# model = 'Europa-tidal1_eta1.0d14_P0.6TW_3.0wt%-NH3_core0.3'
# model = 'Europa-tidal1_eta1.0d14_P0.6TW_3.0wt%-NH3_core0.5'
fig2,(axqq) = plt.subplots(1,1,figsize=(8,6))
axqq.scatter(melting_point_TW,100-melting_point_frac*100,color='#35838D')
axqq.set_ylim(0,100)
axqq.set_xlim(0,1.5)
axqq.set_xlabel('P$_{tide}$ (TW)',fontsize=labelsize)
axqq.set_ylabel('Fraction of tidal heating within ice shell (%)',fontsize=labelsize-2)
rainbow = cm.get_cmap('jet',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
# fig,(ax) = plt.subplots(1,1,figsize=(8,12))
# for i, model in enumerate(model_list):
#     i = i
#     line = open(workpath+model+'_T-profiles.dat').readlines(0)
#     time = line[1].split("at ")[1].split('\n')[0].split('   ')[1:]
#     data = np.loadtxt(workpath+model+'_T-profiles.dat',skiprows=2)
#     mm=30
#     kk = np.hsplit(data, mm+1)
#     uu=pd.DataFrame(kk[mm], columns = ['radius','temperature'])
#     uu.to_csv(path+model+'_Tprofiles_'+time[mm], sep=',', index=False)
    
#     dd = pd.read_csv(path+model+'_Tprofiles_'+time[mm])
        
#     ax.plot(dd.temperature, dd.radius,color=newcolors[i],lw=3)
# ax.set_ylim(0,1561)
# ax.set_xlim(0,2500)
# Pressure = np.linspace(0.21,6.2)
# Tsolidus = 1559.27+42.3 * Pressure
# # Tsolidus2 = 1614.9 + 23.0*Pressure
# depth = Pressure*1e9/3013/1000/1.315
# ax.plot(Tsolidus,1561-depth, c='gray',lw=4,linestyle='dashed')
# # ax.plot(Tsolidus2,1561-depth, c='gray',lw=4,linestyle='dotted')
# # ax.set_title(model,fontsize=labelsize)
# # colorbar
# ax4  = fig.add_axes([0.95,0.12,0.03,0.76])
# norm = mpl.colors.Normalize(vmin=0,vmax=4.55)
# cb1  = mpl.colorbar.ColorbarBase(ax4,cmap=rainbow,norm=norm,orientation='vertical')
# cb1.set_label('Time',fontsize=labelsize)
# cb1.ax.tick_params(axis='y', labelsize=labelsize-2)
# cb1.ax.yaxis.set_label_position('right')
# # ax.set_ylim(200,500)
for aa in [axqq]:
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True,top=True,direction='in',pad=15)
    #aa.set_xlim(xmin,xmax)
    #aa.set_ylim(-zmin,-zmax)
    #aa.set_xlim(0,5)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig2.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figureS7_v1.pdf')
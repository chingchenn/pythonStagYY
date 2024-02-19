#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:59:51 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
temp_profile = 0
plot_Ram_Ftop = 0
variation_P = 1
variation_eta = 1
variation_P_purewater = 0
variation_eta_purewater = 0
labelsize = 30
bwith = 3
ice_core_fraction = 0
tidal_heat = 0
fixP_changeVis_P3 = 0
fixVis_changeP = 0
power_period = 0

### PATH ###
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/StagYY_Works/AGU2023/'

newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']

colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67']


if temp_profile:
    model = 'i10'
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
    end =500
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(12,12))
    # for kk, model in enumerate(model_list):
    ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(end)+'.txt',
                      sep = '\\s+',header = None,names = header_list)    
    ax.plot(ff.vzabs*ff.Tmean,ff.r,color = '#414F67',lw=5)
    ax2.plot(ff.Tmean,ff.r,color = '#414F67',lw=5)
        
    for aa in [ax,ax2]:
        aa.tick_params(labelsize=labelsize)
        for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
        aa.grid()
        aa.set_ylim(0,1)
    ax2.set_xlim(0.0,1)
    ax.set_xlim(0,500)
    
    ax.set_xlabel('T Vz',fontsize = labelsize)
    ax2.set_xlabel('Temperature ',fontsize = labelsize)
    ax.set_ylabel('Depth',fontsize = labelsize)
    # ax.set_title('Snapshot = '+str(end),fontsize = labelsize)
        

    lid_thickness = np.zeros(end)
    new_thickness = np.zeros(end)
    avg_temp = np.zeros(end)
    time = np.linspace(1,end,end)
    for i in range(1,end+1):
        ff = pd.read_csv(modelpath+model+'/datafile/'+model+'_data_'+str(i)+'.txt',
                          sep = '\\s+',header = None,names = header_list)    
        x = np.array(ff.vzabs*ff.Tmean)
        y = np.array(ff.r)
        smooth_d2 = np.gradient(np.gradient(x))
        infls = np.where(np.diff(np.sign(smooth_d2)))[0]
     
        if len(x[infls])==0:
            inf_point = 0
            lid_thickness[i-1]=0
            avg_temp[i-1] = 0
        else:
            inf_point = x[infls][-1]
            index_inf_point = [i for i, kk in enumerate(x) if kk == inf_point][0]
            a = (y[index_inf_point]-y[index_inf_point+1])/(x[index_inf_point]-x[index_inf_point+1])
            b = y[infls][-1]-a*inf_point
            line_x = np.linspace(0,2400)
            line_y = a*line_x + b
            new_thickness[i-1] = line_y[0]
            avg_temp[i-1] = np.average(ff.Tmean[y<y[x==inf_point]])
            
    ax.plot(line_x,line_y,color ='#4198B9',lw = 3 )
    ax.scatter(x[index_inf_point],y[index_inf_point],color = 'orange',s = 200)
    #ax.set_xlim(0,1000)
    ax.axhline(y=line_y[0], color='r', linestyle='--',lw = 2)
    ax2.axhline(y=line_y[0], color='r', linestyle='--',lw = 2)
    ax.axhspan(line_y[0], 1, facecolor='#44b14e',alpha=0.25)
    ax.axhspan(0,line_y[0], facecolor='#32aae2',alpha=0.25)
    ax2.axhspan(line_y[0], 1, facecolor='#44b14e',alpha=0.25)
    ax2.axhspan(0,line_y[0], facecolor='#32aae2',alpha=0.25)
    # fig.savefig(figpath+'temperature_profile.pdf')

if plot_Ram_Ftop:
    a = 1.33
    b = 0.33
    c = 1.70
    
    sigma_a = 0.2
    sigma_b = 0.01
    sigma_c = 0.08
    sigma_F = 0.1
    fig6,(ax) = plt.subplots(1,1,figsize=(12,9))
    labelsize=26
    file = 'plot_Ram_v2013_result_all.dat'
    Ra0_array, f_array, H_array, Deta_array, Tm_array, Ftop_array = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/'+file).T
    # file = 'plot_Ram_v2013_result.dat'
    # Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, eee2 = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/'+file).T
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        H = H_array[kk]
        Deta = np.exp(Deta_array[kk])
        Ftop=Ftop_array[kk]   
        #Ftop = f**2 * Fbot + (1+f+f**2)/3 * H
        gamma=np.log(Deta)
        
    
        Tm = Tm_array[kk] 
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Ftop*(gamma**c)
        x = Ram
        yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
        if H ==0:
            ax.errorbar(x,y,yerr = yerror,fmt='ok',ecolor="#282130") 
        else:
            ax.errorbar(x,y,yerr = yerror,fmt='ok',color='orange',ecolor="#282130") 
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax.plot(x_pred,F_pred,c='#CD5C5C')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax.plot(x_pred,error_plus,c = '#849DAB',lw = 1, linestyle = 'dashed')
    ax.plot(x_pred,error_minus,c = '#849DAB',lw = 1, linestyle = 'dashed')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e5,1e8)
    ax.set_ylim(50,1e3)
    #ax.legend(fontsize=labelsize)
    ax.set_xlabel('Ra$_{eff}$',fontsize = labelsize)
    ax.set_ylabel('$\phi_{top}$',fontsize = labelsize)
    # ax.tick_params(labelsize=labelsize)
    ax.tick_params(labelsize=labelsize,width=1,length=10,right=True, top=True,direction='in', which='both')
    # ax.set_title('normal models with ftop data and ftop inverse',fontsize = 16)
    # ax.grid()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    fig6.savefig(figpath+'Ram_surface heat_flux.pdf')
if variation_P:    
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,2,figsize=(24,20))
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
                ]
    
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e13_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta3e13_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e14_P0.0-2.0TW_3.0%-NH3_k2.6_wt%',]
    
    label_list=['eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint/1e12
        x = Pint
        ax[0].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=4,)
        ax2[0].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=4,)
        ax3[0].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i],lw=4,)
        ax[0].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[0].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[0].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        if len(data.zbot[data.conv==0])>0:
            ax[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint/1e12
        x = Pint
        ax[1].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=4,)
        ax2[1].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=4,)
        ax3[1].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i],lw=4,)
        ax[1].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[1].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[1].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        if len(data.zbot[data.conv==0])>0:
            ax[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax[0].set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2[0].set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3[0].set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax[0].set_ylim(161,0)
    ax2[0].set_ylim(180,280)
    ax3[0].set_ylim(0,40)
    ax[1].set_ylim(161,0)
    ax2[1].set_ylim(180,280)
    ax3[1].set_ylim(0,40)
    ax[0].legend(fontsize = 28)
    ax3[0].set_xlabel('P$_{tide}$',fontsize=labelsize)
    ax3[1].set_xlabel('P$_{tide}$',fontsize=labelsize)
    ax[0].set_title('2D scaling parameters',fontsize=labelsize)
    ax[1].set_title('3D scaling parameters',fontsize=labelsize)
    for aa in [ax[0],ax2[0],ax3[0],ax[1],ax2[1],ax3[1]]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(0.0,2)
        aa.grid()
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)

if variation_P_purewater:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,2,figsize=(24,20))
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
                ]
    
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e13_P0.0-2.0TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta3e13_P0.0-2.0TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e14_P0.0-2.0TW_0.0%-NH3_k2.6_wt%',]
    
    label_list=['eta = 1e13','eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint/1e12
        x = Pint
        ax[0].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=4,)
        ax2[0].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=4,)
        ax3[0].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i],lw=4,)
        ax[0].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[0].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[0].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        if len(data.zbot[data.conv==0])>0 and len(data.zbot[data.conv==1])>0:
            ax[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        Pint=data.Pint/1e12
        x = Pint
        ax[1].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i],lw=4,)
        ax2[1].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],lw=4,)
        ax3[1].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i],lw=4,)
        ax[1].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[1].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[1].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        if len(data.zbot[data.conv==0])>0:
            ax[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax[0].set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2[0].set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3[0].set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax[0].set_ylim(161,0)
    ax2[0].set_ylim(180,280)
    ax3[0].set_ylim(0,40)
    ax[1].set_ylim(161,0)
    ax2[1].set_ylim(180,280)
    ax3[1].set_ylim(0,40)
    ax[0].legend(fontsize = 28)
    ax3[0].set_xlabel('P$_{tide}$',fontsize=labelsize)
    ax3[1].set_xlabel('P$_{tide}$',fontsize=labelsize)
    ax[0].set_title('2D scaling parameters',fontsize=labelsize)
    ax[1].set_title('3D scaling parameters',fontsize=labelsize)
    for aa in [ax[0],ax2[0],ax3[0],ax[1],ax2[1],ax3[1]]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(0.0,2)
        aa.grid()
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
if variation_eta:    
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,2,figsize=(24,20))
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
                ]
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.3TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.0TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.2TW_3.0%-NH3_k2.6_wt%',]
                # 'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.0TW_3.0%-NH3_k2.6_wt%']
    
    label_list=['P = 0.0 Tw','P = 0.3 Tw','P = 0.6 Tw','P = 1.0 Tw','P = 1.2 Tw',]
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x=data.etaref
        ax[0].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i], lw=4,)
        ax2[0].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i], lw=4,)
        ax3[0].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i], lw=4,)
        ax[0].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[0].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[0].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        if len(data.zbot[data.conv==0])>0:
            ax[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x=data.etaref
        ax[1].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i], lw=4,)
        ax2[1].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i], lw=4,)
        ax3[1].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i], lw=4,)
        ax[1].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[1].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[1].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        if len(data.zbot[data.conv==0])>0:
            ax[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax[0].set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2[0].set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3[0].set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax[0].set_ylim(161,0)
    ax2[0].set_ylim(180,280)
    ax3[0].set_ylim(0,40)
    ax[1].set_ylim(161,0)
    ax2[1].set_ylim(180,280)
    ax3[1].set_ylim(0,40)
    ax[0].legend(fontsize = 25)
    # ax[0].legend(fontsize = 18)
    ax[0].set_title('2D scaling parameters',fontsize=labelsize)
    ax[1].set_title('3D scaling parameters',fontsize=labelsize)
    ax3[0].set_xlabel('log(reference viscosity)',fontsize=labelsize)
    ax3[1].set_xlabel('log(reference viscosity)',fontsize=labelsize)
    for aa in [ax[0],ax2[0],ax3[0],ax[1],ax2[1],ax3[1]]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.grid()
        aa.set_xlim(1e12,1e15)
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    # fig.savefig(figpath+'scaling_eta_NH3.pdf')
if variation_eta_purewater:
    header_list = ['Pint','etaref','pCLA','conv','P','zbot','%vol',
                   'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
    fig,(ax,ax2,ax3) = plt.subplots(3,2,figsize=(24,20))
    model_list=['xtest-Europa_eta1e14_P0.0-0.3TW_3.0%-NH3_k2.6_wt%',
                ]
    model_list=['Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.0TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.3TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P0.6TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.0TW_0.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.2TW_0.0%-NH3_k2.6_wt%',]
                # 'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e12-1e15_P1.0TW_3.0%-NH3_k2.6_wt%']
    
    label_list=['P = 0.0 Tw','P = 0.3 Tw','P = 0.6 Tw','P = 1.0 Tw','P = 1.2 Tw',]
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x=data.etaref
        ax[0].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i], lw=4,)
        ax2[0].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i], lw=4,)
        ax3[0].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i], lw=4,)
        ax[0].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[0].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[0].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        if len(data.zbot[data.conv==0])>0:
            ax[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[0].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    for i, model in enumerate(model_list):
        workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/old_scaling_thermal_evolution/'
        data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x=data.etaref
        ax[1].plot(x[data.conv==1], data.zbot[data.conv==1],color=colors[i],label = label_list[i], lw=4,)
        ax2[1].plot(x[data.conv==1],data.Tm[data.conv==1],color=colors[i], lw=4,)
        ax3[1].plot(x[data.conv==1], data.dlid[data.conv==1],color=colors[i], lw=4,)
        ax[1].plot(x[data.conv==0], data.zbot[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax2[1].plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3[1].plot(x[data.conv==0], data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        if len(data.zbot[data.conv==0])>0:
            ax[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                       ymax=data.zbot[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax2[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.Tm[data.conv==0].iloc[0], 
                       ymax=data.Tm[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
            ax3[1].vlines(x=x[data.conv==0].iloc[0], ymin=data.dlid[data.conv==0].iloc[0], 
                       ymax=data.dlid[data.conv==1].iloc[-1], colors=colors[i], ls='--', lw=3,)
    ax[0].set_ylabel('ice layer thickness (km)',fontsize = labelsize)
    ax2[0].set_ylabel('interior temperature (K)',fontsize = labelsize-5)
    ax3[0].set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
    ax[0].set_ylim(161,0)
    ax2[0].set_ylim(180,280)
    ax3[0].set_ylim(0,40)
    ax[1].set_ylim(161,0)
    ax2[1].set_ylim(180,280)
    ax3[1].set_ylim(0,40)
    ax[0].set_title('2D scaling parameters',fontsize=labelsize)
    ax[1].set_title('3D scaling parameters',fontsize=labelsize)
    ax3[0].set_xlabel('log(reference viscosity)',fontsize=labelsize)
    ax3[1].set_xlabel('log(reference viscosity)',fontsize=labelsize)
    for aa in [ax[0],ax2[0],ax3[0],ax[1],ax2[1],ax3[1]]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.grid()
        aa.set_xlim(1e12,1e15)
        aa.set_xscale('log')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    # fig.savefig(figpath+'scaling_eta_water.pdf')
if ice_core_fraction:
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']
    model_list=['Europa-tidal1_period0.0Gyr_m0.0_core1.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.2_core0.8_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.5_core0.5_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.7_core0.3_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%']

    model_list=['Europa-tidal1_period0.0Gyr_m0.0_core1.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.25_core0.75_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.5_core0.5_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.75_core0.25_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',]
    model_list=['Europa-tidal1_period0.0Gyr_m0.0_core1.0_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.25_core0.75_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.5_core0.5_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m0.75_core0.25_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                'Europa-tidal1_period0.0Gyr_m1.0_core0.0_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['Europa-tidal1_period0.0Gyr_m0.75_core0.25_eta3e13_P0.6TW_3.0%-NH3_k2.6_wt%',]
    
    label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0 Tw','P = 1 Tw','P = 2 Tw','P = 3 Tw',]
    label_list = ['0.125 Gyr','0.25 Gyr','0.5 Gyr','1.0 Gyr']
    label_list = ['0 %','25 %','50 %','75 %','100 %',]
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(20,18))
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        for kk in range(3,len(data)-1):
            if data.conv[kk-1]*data.conv[kk]==1 and data.conv[kk]*data.conv[kk+1]==0:
                print(data.time_Gyr[kk],'conv-->cond')
            if data.conv[kk-1]*data.conv[kk]==0 and data.conv[kk]*data.conv[kk+1]==1:
                print(data.time_Gyr[kk],'cond-->conv')
            # if data.conv[kk]==1:
            #     ax.plot(x[kk],data.zbot[kk],color=colors[i],lw=3)
            #     ax2.plot(x[kk],data.Tm[kk],color=colors[i],lw=3,label=label_list[i])
            #     ax3.plot(x[kk],data.dlid[kk],color=colors[i],lw=3)
            # if data.conv[kk]==0:
            #     ax.plot(x[kk], data.zbot[kk],color=colors[i],linestyle='dashed',lw=3)
            #     ax2.plot(x[kk],data.Tm[kk],color=colors[i],linestyle='dashed',lw=3)
            #     ax3.plot(x[kk], data.dlid[kk],color=colors[i],linestyle='dashed',lw=3)
        ax.scatter(x[data.conv==1],data.zbot[data.conv==1],color=colors[i],s=20)
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i],s=20,label=label_list[i])
        ax3.scatter(x[data.conv==1],data.dlid[data.conv==1],color=colors[i],s=20)
        ax.plot(x[data.conv==0],data.zbot[data.conv==0],color='gray')
        ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color='gray')
        ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color='gray')
        
        
        # ax.plot(x,data.Pint,color=colors[i],lw=3)
        
        ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
        ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
        ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
        ax.set_ylim(161,0)
        ax2.set_ylim(180,280)
        ax3.set_ylim(0,70)
        ax2.legend(fontsize=labelsize-10)
        for aa in [ax,ax2,ax3]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
            ax3.set_xlabel('Time',fontsize=labelsize)
            aa.set_xlim(0,4.55)
            aa.grid()
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
if tidal_heat:
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
    model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',]
    model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',]
    
    model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',]
    
    model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',]

    # model_list = ['Europa-tidal5_period0.125Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',]
    # # model_list = ['Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
    # #                'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%']
    model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_dcay40%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_dcay70%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%']
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.6TW_3.0%-NH3_k2.6_wt%_Ea40',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_Ea40',
    #               ]
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%_Ea40',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%_Ea40',
    #               ]
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.1TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',]
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',]
    
    # model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%_Ea40',
    #               'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',]
    colors=['#2F4F4F','#4682B4','#CD5C5C','#97795D',]
    fig2,(ax4) = plt.subplots(1,1,figsize=(18,8))
    label_list=['eta = 3.2e13','eta = 1e14','eta = 3.2e14']
    label_list=['P = 0.3 Tw','P = 0.6 Tw','P = 2 Tw','P = 3 Tw',]
    label_list = ['0.125 Gyr','0.25 Gyr','0.5 Gyr','1.0 Gyr']
    label_list = ['period of 1.00 Gyr with 20% decay','period of 1.00 Gyr with 40% decay',
                  'period of 1.00 Gyr with 70% decay',
                  'period of 0.25 Gyr with 20% decay','0.5 Gyr','1.0 Gyr']
    
    
    
    
    
    
    # from scipy.datasets import electrocardiogram
    from scipy.signal import find_peaks
    fig,(ax,ax2,ax3) = plt.subplots(3,1,figsize=(20,18))
    for i, model in enumerate(model_list):
        # i = 3
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        to_cond_time_list = []
        end_cond_time_list = []
        for kk in range(3,len(data)-1):
            if data.conv[kk-1]*data.conv[kk]==1 and data.conv[kk]*data.conv[kk+1]==0:
                # print(data.time_Gyr[kk],'conv-->cond')
                to_cond_time_list.append(data.time_Gyr[kk])
            if data.conv[kk-1]*data.conv[kk]==0 and data.conv[kk]*data.conv[kk+1]==1:
                # print(data.time_Gyr[kk],'cond-->conv')
                end_cond_time_list.append(data.time_Gyr[kk])
        time_interval = []
        for mm in range(len(to_cond_time_list)):
            # print(end_cond_time_list[mm]-to_cond_time_list[mm])
            time_interval.append(end_cond_time_list[mm]-to_cond_time_list[mm])
        print(np.average(time_interval))
        
        
        
        # ax.scatter(x[data.conv==1], data.zbot[data.conv==1])
        ax.scatter(x[data.conv==1],data.zbot[data.conv==1],color=colors[i])
        ax2.scatter(x[data.conv==1],data.Tm[data.conv==1],color=colors[i])#,label=label_list[i])
        ax3.scatter(x[data.conv==1],data.dlid[data.conv==1],color=colors[i],label=label_list[i])
        
        ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=colors[i],linestyle='dashed',label=label_list[i])
        ax2.plot(x[data.conv==0],data.Tm[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        ax3.plot(x[data.conv==0],data.dlid[data.conv==0],color=colors[i],linestyle='dashed',lw=3)
        
        # x = electrocardiogram()[2000:4000]
        peaks, _ = find_peaks(data.zbot[data.conv==1], height=0)
        print(data.zbot[data.conv==1].iloc[peaks])
        # plt.plot(x)
        # plt.plot(peaks, x[peaks], "x")
        # ax2.plot(x,data.Pint,color=colors[i],lw=3)
        ax4.plot(x,data.Pint,color=colors[i],lw=3,label=label_list[i])

        
        ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
        ax2.set_ylabel('interior temperature (K)',fontsize = labelsize-5)
        ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
        ax.set_ylim(161,0)
        ax3.set_ylim(50,0)
        ax3.set_xlabel('Time',fontsize=labelsize)
        ax4.set_xlabel('Time',fontsize=labelsize)
        ax4.set_ylabel('Power (TW)',fontsize=labelsize)
        ax4.set_ylim(0,4)
        ax4.legend(fontsize=labelsize)
        ax3.legend(fontsize=labelsize-5)
        for aa in [ax,ax2,ax3,ax4]:
            aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
            aa.set_xlim(0,4.55)
            # aa.grid()
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
        # fig2.savefig(figpath+'Tidal heating.pdf')
        # fig3,(ax1) = plt.subplots(1,1,figsize=(18,10))
        # ax1.bar(x[data.conv==1].iloc[peaks],data.zbot[data.conv==1].iloc[peaks],width=0.17,color='seagreen',label='peridotite')
    
        # for aa in [ax1]:
        #     aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        #     aa.set_xlim(0,4.55)
        #     # aa.grid()
        #     for axis in ['top','bottom','left','right']:
        #         aa.spines[axis].set_linewidth(bwith)

colors = ['#2F4F4F','#4682B4','#CD5C5C','#97795D',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

if fixP_changeVis_P3:
    colors = ['#2F4F4F','#4682B4','#CD5C5C','#97795D',
                  '#6B0D47','#282130','#7E9680','#24788F',
                  '#849DAB','#EA5E51','#35838D','#4198B9',
                  '#414F67','#97795D','#6B0D47',]
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/thermal_evolution/'
    model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta6e12_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta3e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta6e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e14_P0.3TW_3.0%-NH3_k2.6_wt%',]
    header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
                   'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']
    
    label_list=['eta = 6e12','eta = 1e13','eta = 3e13','eta = 6e13','eta = 1e14','eta = 3.2e14']
    fig,(ax,ax3) = plt.subplots(2,1,figsize=(15,10),gridspec_kw={'height_ratios':[1,0.7]})
    # ax.set_title('Europa-tidal5_period1.0Gyr_m1.0_core0.0_etaXXX_P0.3TW_3.0%-NH3_k2.6_wt%',fontsize=labelsize)
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        
        zbot_cond =ma.array(data.zbot, mask = mask_cond)
        zbot_conv=ma.array(data.zbot, mask = mask_conv)
        dlid_cond = ma.array(data.dlid, mask = mask_cond)
        dlid_conv = ma.array(data.dlid, mask = mask_conv)
        
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.plot(x,dlid_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.axhline(y=15.5, color='gray', linestyle='--',lw = 2)
        # print(model,np.min(data.zbot[(data.conv==1)*(data.time_Gyr>1)]))
        # print(model,np.max(data.dlid[(data.conv==1)*(data.time_Gyr>1)]))
        
        ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
        # ax3.plot(x,dlid_cond,color=colors[i],linestyle='dashed')
        ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize-5)
        ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
        ax.set_ylim(161,0)
        ax3.set_ylim(25,0)
        ax3.set_xlabel('time',fontsize=labelsize-3)
        # ax3.legend(fontsize=labelsize-10)
        for aa in [ax,ax3]:
            aa.tick_params(labelsize=labelsize-10,width=3,length=10,right=True, top=True,direction='in',pad=5)
            aa.set_xlim(0,4.55)
            # aa.grid()
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
    # fig.savefig(figpath+'thermal_model_varing_viscosity.pdf')

if fixVis_changeP:
    colors=['#282130','#849DAB','#35838D','#6B0D47','#97795D','#414F67']
    model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.1TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.3TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P1.0TW_3.0%-NH3_k2.6_wt%',]
                  #'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P1.2TW_3.0%-NH3_k2.6_wt%',]
    
    label_list=['P = 0.1 Tw','P = 0.3 Tw','P = 0.6 Tw','P = 1.0 Tw','P = 1.2 Tw',]
    fig,(ax,ax3) = plt.subplots(2,1,figsize=(15,10),gridspec_kw={'height_ratios':[1,0.7]})
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        
        zbot_cond =ma.array(data.zbot, mask = mask_cond)
        zbot_conv=ma.array(data.zbot, mask = mask_conv)
        dlid_cond = ma.array(data.dlid, mask = mask_cond)
        dlid_conv = ma.array(data.dlid, mask = mask_conv)
        
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.plot(x,dlid_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.axhline(y=6.3, color='gray', linestyle='--',lw = 2)
        
        ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
        ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize-5)
        ax3.set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
        ax.set_ylim(161,0)
        ax3.set_ylim(25,0)
        ax3.set_xlabel('time',fontsize=labelsize)
        ax3.legend(fontsize=labelsize-10)
        for aa in [ax,ax3]:
            aa.tick_params(labelsize=labelsize-10,width=3,length=10,right=True, top=True,direction='in',pad=5)
            aa.set_xlim(0,4.55)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
    fig.savefig(figpath+'thermal_model_varing_power.pdf')

if power_period:
    colors = ['#2F4F4F','#97795D','#4682B4','#CD5C5C',
                  '#6B0D47','#282130','#7E9680','#24788F',
                  '#849DAB','#EA5E51','#35838D','#4198B9',
                  '#414F67','#97795D','#6B0D47',]
    model_list = ['Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.25Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',
                  'Europa-tidal5_period0.5Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_dcay40%',
                  #'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_dcay40%',
                  'Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%_dcay70%',
                  ]
    
    label_list=['period of 1.00 Gyr with 20% decay','period of 0.25 Gyr with 20% decay',
                'period of 0.5 Gyr with 40% decay',
                #'period of 1.00 Gyr with 40% decay', 
                'period of 1.00 Gyr with 70% decay',]
    fig2,(ax4) = plt.subplots(1,1,figsize=(20,6))
    # ax.set_title('Europa-tidal5_period1.0Gyr_m1.0_core0.0_eta1e13_P0.6TW_3.0%-NH3_k2.6_wt%',fontsize=labelsize)
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        ax4.plot(x,data.Pint,color=colors[i],lw=4,label=label_list[i])
        ax4.set_ylabel('power (TW)',fontsize=labelsize)
        ax4.set_ylim(0,3)
        # ax4.legend(fontsize=labelsize)
        ax4.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        ax4.set_xlim(0,4.55)
        for axis in ['top','bottom','left','right']:
            ax4.spines[axis].set_linewidth(bwith)
        ax4.set_xlabel('time',fontsize=labelsize)
    # fig2.savefig(figpath+'period_power.pdf')
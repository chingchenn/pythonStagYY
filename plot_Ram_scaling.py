#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 23:04:55 2023

@author: chingchen
"""

import numpy as np
import matplotlib.pyplot as plt


workpath = '/Users/chingchen/Desktop/StagYY_Works/'

fig5 = 0 ## normal models top
fig6 = 1 ## H models top
fig7 = 0 ## normal models bottom
fig8 = 0 ## H models bottom
heat_bot = 0
heat_top = 1
if heat_bot:
    a = 1.5448
    b = 0.3023
    c = 1.4697
    
    sigma_a = 0.3182
    sigma_b = 0.01213
    sigma_c = 0.1195
    sigma_F = 0.02
if heat_top:   
    a = 1.33
    b = 0.33
    c = 1.70
    
    sigma_a = 0.2
    sigma_b = 0.01
    sigma_c = 0.08
    sigma_F = 0.1
    
    # a = 1.46
    # b = 0.27
    # c = 1.21
    
    # sigma_a = 0.2
    # sigma_b = 0.01
    # sigma_c = 0.01
    # sigma_F = 0.2
   
if fig5:
    fig5,(ax) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    path = '/Users/chingchen/Desktop/StagYY_Works/nlg_inverse_3v3p_solve_Tm/'
    #file = 'scalinglaw_rawdata_ftop.txt'
    #file = 'scalinglaw_inversion_jiching_top_1e6_read'
    #file = 'normal_inversion_read.dat'
    file = 'normal_inversion_1115.dat'
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, eee2 = np.loadtxt(path+file).T
    
    #Ra0_array, f_array, H_array, Deta_array, Tm_array, Ftop_array = np.loadtxt(path+file).T
    
    # npt=1
    
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        Deta = np.exp(Deta_array[kk])
        Ftop=Ftop_array[kk]   
        gamma=np.log(Deta)
        
    
        Tm = Tm_array[kk] 
        
        
        
        
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Ftop*(gamma**c)
        x = Ram
        #print(Ram)
        
        yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
        ax.errorbar(x,y,yerr = yerror,fmt='o',c='#AE6378',ecolor="darkblue") 
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.legend(fontsize=labelsize)
    ax.set_xlabel('Ram',fontsize = labelsize)
    ax.set_ylabel('F_top',fontsize = labelsize)
    ax.tick_params(labelsize=labelsize)
    #ax.set_xticks(np.logspace(3,5,3))
    # ax.set_xlim(1e5,1e9)
    #ax.set_ylim(20,1200)
    ax.set_title('normal models with ftop data and ftop inverse',fontsize = 16)
    ax.grid()
    ###------------------------------------------
    # predicting 
    a=a+0.2
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')

if fig6:
    fig6,(ax) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    file = 'plot_Ram_v2013_result_all.dat'
    Ra0_array, f_array, H_array, Deta_array, Tm_array, Ftop_array = np.loadtxt(workpath+file).T
    file = 'plot_Ram_v2013_result.dat'
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, eee2 = np.loadtxt(workpath+file).T
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
        ax.errorbar(x,y,yerr = yerror,fmt='ok',ecolor="r") 
        
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e5,1e9)
    ax.set_ylim(1e2,1e3)
    #ax.legend(fontsize=labelsize)
    ax.set_xlabel('Ram',fontsize = labelsize)
    ax.set_ylabel('F_top',fontsize = labelsize)
    ax.tick_params(labelsize=labelsize)
    ax.set_title('normal models with ftop data and ftop inverse',fontsize = 16)
    ax.grid()
####################################################################################################
if fig7:
    fig7,(ax7) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Fbot_array, eee2 = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/scalinglaw_rawdata_fbot.txt').T
    
    
    
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        Deta = np.exp(Deta_array[kk])
        Fbot=Fbot_array[kk]   
        gamma=np.log(Deta)
        Ramed=Ra0*np.exp(0.5*gamma)
    
        Tm = Tm_array[kk] 
        #Ram = Ra0*np.exp(gamma*(Tm-0.5)) 
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Fbot*f*(gamma**c)
        x = Ram
        
        
        yerror = sigma_F/Fbot* y + sigma_c * y * np.log(gamma)
        ax7.errorbar(x,y,yerr = yerror,fmt='og',ecolor="r") 
        print(yerror)
        
        
    ax7.set_yscale('log')
    ax7.set_xscale('log')
    #ax.legend(fontsize=labelsize)
    ax7.set_xlabel('Ram',fontsize = labelsize)
    ax7.set_ylabel('F_bot',fontsize = labelsize)
    ax7.tick_params(labelsize=labelsize)
    #ax.set_xticks(np.logspace(3,5,3))
    ax7.set_xlim(1e5,1e9)
    ax7.set_ylim(40,1200)
    ax7.set_title('normal models with fbot data and fbot inverse',fontsize = 16)
    ax7.grid()
    ###------------------------------------------
    ## predicting 
    
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax7.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax7.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax7.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')
    

if fig8:
    fig8,(ax8) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Fbotp_array, H_array = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/scalinglaw_H_rawdata_fbot.txt').T
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, H_array = np.loadtxt(workpath+'scalinglaw_H_resolution_ftop.txt').T

    npt=1
    
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        H = H_array[kk]
        Deta = np.exp(Deta_array[kk])
        Fbot=Fbotp_array[kk]   
        Ftop = f**2 * Fbot + (1+f+f**2)/3 * H
        gamma=np.log(Deta)
        Ramed=Ra0*np.exp(0.5*gamma)
    
        Tm = Tm_array[kk] 
        #Ram = Ra0*np.exp(gamma*(Tm-0.5)) 
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Ftop/f*(gamma**c)
        x = Ram
        
        
        yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
        ax8.errorbar(x,y,yerr = yerror,fmt='og',ecolor="r") 
        print(format(x,'.1e'))
        
        
    ax8.set_yscale('log')
    ax8.set_xscale('log')
    #ax.legend(fontsize=labelsize)
    ax8.set_xlabel('Ram',fontsize = labelsize)
    ax8.set_ylabel('F_bot',fontsize = labelsize)
    ax8.tick_params(labelsize=labelsize)
    #ax.set_xticks(np.logspace(3,5,3))
    ax8.set_xlim(1e5,1e9)
    ax8.set_ylim(40,1200)
    ax8.set_title('internal heating models with fbot data convert to ftop and fbot inverse',fontsize = 16)
    ax8.grid()
    ###------------------------------------------
    ## predicting 
    
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax8.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax8.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax8.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')



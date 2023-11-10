#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 23:04:55 2023

@author: chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

fig5 = 1
fig6 = 0
fig7 = 0
fig8 = 0
if fig5:
    fig5,(ax) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    path = '/Users/chingchen/Desktop/StagYY_Works/nlg_inverse_3v3p_solve_Tm/'
    file = 'scalinglaw_rawdata_ftop.txt'
    file = 'scalinglaw_inversion_jiching_top_1e6_read'
    #file = 'normal_inversion_read.dat'
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, eee2 = np.loadtxt(path+file).T
    
    heat_bot = 1
    heat_top = 0
    if heat_bot:
        a = 1.5448
        b = 0.3023
        c = 1.4697
        
        sigma_a = 0.3182
        sigma_b = 0.01213
        sigma_c = 0.1195
        sigma_F = 0.02
    if heat_top:
        a = 1.5777
        b = 0.2878
        c = 1.3904
        
        sigma_a = 0.3361
        sigma_b = 0.01388
        sigma_c = 0.1287
        sigma_F = 0.02
    
    npt=1
    
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        Deta = np.exp(Deta_array[kk])
        Ftop=Ftop_array[kk]   
        gamma=np.log(Deta)
        Ramed=Ra0*np.exp(0.5*gamma)
    
        Tm = Tm_array[kk] 
        #Ram = Ra0*np.exp(gamma*(Tm-0.5)) 
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Ftop/f*(gamma**c)
        x = Ram
        
        
        yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
        ax.errorbar(x,y,yerr = yerror,fmt='og',ecolor="r") 
        print(yerror)
        
        
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.legend(fontsize=labelsize)
    ax.set_xlabel('Ram',fontsize = labelsize)
    ax.set_ylabel('F_bot',fontsize = labelsize)
    ax.tick_params(labelsize=labelsize)
    #ax.set_xticks(np.logspace(3,5,3))
    ax.set_xlim(1e5,1e9)
    ax.set_ylim(40,1200)
    ax.set_title('normal models with ftop data and ftop inverse',fontsize = 16)
    ax.grid()
    ###------------------------------------------
    ## predicting 
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')

if fig6:
    fig6,(ax2) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, H_array = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/scalinglaw_H_rawdata_ftop.txt').T
    
    npt=1
    
    Ra_model = np.logspace(4,10)
    for kk in range(len(Ra0_array)):
        Ra0 = Ra0_array[kk]
        f = f_array[kk]
        H = H_array[kk]
        Deta = np.exp(Deta_array[kk])
        Ftop=Ftop_array[kk]   
        #Ftop = f**2 * Fbot + (1+f+f**2)/3 * H
        gamma=np.log(Deta)
        Ramed=Ra0*np.exp(0.5*gamma)
    
        Tm = Tm_array[kk] 
        #Ram = Ra0*np.exp(gamma*(Tm-0.5)) 
        Ram = Ra0*np.exp(gamma*Tm) 
        y = Ftop/f*(gamma**c)
        x = Ram
        
        
        yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
        ax2.errorbar(x,y,yerr = yerror,fmt='og',ecolor="r") 
        #print(yerror)
        
        
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    #ax.legend(fontsize=labelsize)
    ax2.set_xlabel('Ram',fontsize = labelsize)
    ax2.set_ylabel('F_bot',fontsize = labelsize)
    ax2.tick_params(labelsize=labelsize)
    #ax.set_xticks(np.logspace(3,5,3))
    ax2.set_xlim(1e5,1e9)
    ax2.set_ylim(40,1200)
    ax2.set_title('internal heating models with ftop data and ftop inverse',fontsize = 16)
    
    ax2.grid()
    ###------------------------------------------
    ## predicting 
    
    
    F_pred = a * Ra_model **b
    x_pred = Ra_model
    ax2.plot(x_pred,F_pred,c='b')
    error_plus = (a+sigma_a) * Ra_model **(b+sigma_b)
    error_minus = (a-sigma_a) * Ra_model **(b-sigma_b)
    ax2.plot(x_pred,error_plus,c = 'k',lw = 1, linestyle = 'dashed')
    ax2.plot(x_pred,error_minus,c = 'k',lw = 1, linestyle = 'dashed')
    


####################################################################################################
if fig7:
    fig7,(ax7) = plt.subplots(1,1,figsize=(12,9))
    labelsize=16
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Fbot_array, eee2 = np.loadtxt('/Users/chingchen/Desktop/StagYY_Works/scalinglaw_rawdata_fbot.txt').T
    
    heat_bot = 1
    heat_top = 0
    if heat_bot:
        a = 1.5448
        b = 0.3023
        c = 1.4697
        
        sigma_a = 0.3182
        sigma_b = 0.01213
        sigma_c = 0.1195
        sigma_F = 0.02
    if heat_top:
        a = 1.5653
        b = 0.3002
        c = 1.4625
        
        sigma_a = 0.3116
        sigma_b = 0.01218
        sigma_c = 0.1188
        sigma_F = 0.02
    
    
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
        #print(yerror)
        
        
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



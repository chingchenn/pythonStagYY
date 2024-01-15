#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 23:04:55 2023

@author: chingchen
"""

import numpy as np
import matplotlib.pyplot as plt


workpath = '/Users/chingchen/Desktop/StagYY_Works/'

a = 1.33
b = 0.33
c = 1.70

sigma_a = 0.2
sigma_b = 0.01
sigma_c = 0.08
sigma_F = 0.1

a = 1.934
b = 0.339
c = 1.87

sigma_a = 0.6
sigma_b = 0.01
sigma_c = 0.1
sigma_F = 0.01

fig6,(ax) = plt.subplots(1,1,figsize=(12,9))
labelsize=16
file = 'plot_Ram.dat'
Ra0_array, f_array, H_array, Deta_array, Tm_array, Ftop_array = np.loadtxt(workpath+file).T
Ra_model = np.logspace(4,10)
for kk in range(len(Ra0_array)):
    Ra0 = Ra0_array[kk]
    f = f_array[kk]
    H = H_array[kk]
    Deta = np.exp(Deta_array[kk])
    Ftop=Ftop_array[kk]   
    gamma=np.log(Deta)
    
    Tm = Tm_array[kk] 
    Ram = Ra0*np.exp(gamma*Tm) 
    y = Ftop*(gamma**c)
    x = Ram
    
    yerror = sigma_F/Ftop* y + sigma_c * y * np.log(gamma)
    if H_array[kk]>0:
        # ax.errorbar(x,y,yerr = yerror,fmt='ok',ecolor="r",label='mix heat' if kk==12 else '') 
        ax.errorbar( x,y,yerr = yerror, fmt='o', capsize=5,color = 'orange',ecolor='#849DAB',label='mix heat' if kk==12 else '')
    else:
        # ax.errorbar(x,y,yerr = yerror,fmt='og',ecolor="r",label='basal heat' if kk==0 else '') 
        ax.errorbar(x,y,yerr = yerror, fmt='o', capsize=5,color = 'k',ecolor='#849DAB',label='basal heat' if kk==0 else '')
    
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
ax.set_ylim(1e1,5e3)
#ax.legend(fontsize=labelsize)
ax.set_xlabel('Ram',fontsize = labelsize)
ax.set_ylabel('F_top',fontsize = labelsize)
ax.tick_params(labelsize=labelsize)
ax.set_title('normal models with ftop data and ftop inverse',fontsize = 16)
ax.grid()
ax.legend(fontsize = labelsize)
####################################################################################################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:47:20 2023

@author: chingchen
"""
import numpy as np
import matplotlib.pyplot as plt


normal   = 0
internal = 1
fig5,(ax,ax2) = plt.subplots(1,2,figsize=(18,9))
fig6,(ax3) = plt.subplots(1,1,figsize=(12,9))
def solve_temperature(a1,a2,c,d,icalc,itemax,nite,eps,f,cgeom,H,gamma,RT,Rk,Ra0,T_TDV,Tm,Ram,km):
    T_H = (a1-a2*f)*(cgeom*H)**c / Ram**d
    g =Tm-T_TDV -T_H
    test = abs(g)
    ite=1
    while (test>eps) and (ite < itemax):  
        deriv = 1 + T_H*d*(gamma + RT*km)
        ordin = g - deriv*Tm
        Tm = -ordin/deriv
        km = 1.0/(1.0 + RT*Tm)
        Ram = Ra0*np.exp(gamma*Tm)/km
        T_H = (a1 - a2*f) * (cgeom*H)**c / Ram**d
        g = Tm - T_TDV - T_H
        test=abs(g)
        print(ite,Tm,Ram,test)
        ite+=1
        
    if icalc ==1:
        nite = ite-1
    return ite,Tm,Ram,test
#-----------------------------------------------------------------------------------|
#  Calculates internal temperature in a system animated by stagnant-lid convection. |
#  Takes into account mixed-heating and temperature-dependent thermal conductivity, |
#  if requested.                                                                    |
#                                                                                   |
#  This version allows either a single calculation (icalc = 1), or calculation of   |
#  temperaturem as a function of one entry parameter: surface Ra (Rasurf), rate of  |
#  heating, H, Thermal viscosity ratio, Deta, or top-to-bottom conductivity ratio,  |
#  Rk. See input files for more details.                                            |
#                                                                                   |
#  Compile with fortran 95:                                                         |  
#              gfortran -o solve_Tm_TDV-H solve_Tm_TDV-H_mod.f95 solve_Tm_TDV-H.f95 |
#                                                                                   |
#  Frederic Deschamps, 2020                                                         |
#  Ji-Ching Chen, 2023                                                              |
#-----------------------------------------------------------------------------------|

path = '/Users/chingchen/Desktop/StagYY_Works/nlg_inverse_3v3p_solve_Tm/'

if normal:
    file = 'normal_for_Tm.dat'
    Ra0_array, Deta_array, f_array, Tm_array, eee1, Ftop_array, eee2 = np.loadtxt(path+file).T
    H = 0.0
    Rk = 1
    
if internal:
    file = 'internal_for_Tm_jiching.dat'
    Ra0_array, Deta_array, f_array, H_array,Tm_array, Ftop_array, Fbot_array= np.loadtxt(path+file).T
    Rk = 1
npt=1
    
#print('Parameters a, b and e of T_TDV scaling, a*Rk^e/gamma/f^b :')
a,b,e = 1.23,1.5,0.2
#print('Uncertainties:')
sigma_a,sigma_b,sigma_e=0.05,0.02,0.02
#print('Parameters a1, a2, c and d of T_H scaling (a1-a2*f)*(cgeom*H)^c/Ram^d for Urey < 1:')
a1_Fpos,a2_Fpos,c_Fpos,d_Fpos=6.0359,5.12078,1.00,0.250
#print('Uncertainties:')
sigma_a1_Fpos,sigma_a2_Fpos,sigma_c_Fpos,sigma_d_Fpos=0.3929,0.5261,0.001,0.001
#print('Parameters a1, a2, c and d of T_H scaling (a1-a2*f)*(cgeom*H)^c/Ram^d for Urey > 1:')
a1_Fneg,a2_Fneg,c_Fneg,d_Fneg=5.36,3.00,1.72,0.333
#print('Uncertainties:')
sigma_a1_Fneg,sigma_a2_Fneg,sigma_c_Fneg,sigma_d_Fneg=0.15,0.15,0.0,0.0
#print( 'Parameters ax, bx, cx, dx, ex of top heat flux scaling, Ftop = ax*Ram^bx*f^dx/gamma^cx/Rk^ex for Urey < 1:')
aphi_Fpos,bphi_Fpos,cphi_Fpos,dphi_Fpos,ephi_Fpos=1.934,0.339,1.87,0.0,0.82
#print('Uncertainties:')
sigma_aphi_Fpos,sigma_bphi_Fpos,sigma_cphi_Fpos,sigma_dphi_Fpos,sigma_ephi_Fpos=0.06,0.004,0.03,0.0,0.01
#print('Parameters ax, bx, cx, dx, ex of top heat flux scaling, Ftop = ax*Ram^bx*f^dx/gamma^cx/Rk^ex for Urey > 1:')
aphi_Fneg,bphi_Fneg,cphi_Fneg,dphi_Fneg,ephi_Fneg=1.57,0.270,1.21,0.0,0.82
#print('Uncertainties:')
sigma_aphi_Fneg,sigma_bphi_Fneg,sigma_cphi_Fneg,sigma_dphi_Fneg,sigma_ephi_Fneg=0.06,0.004,0.03,0.0,0.01
#print('Threshold value, maximum number of iterations')
eps,itemax=1.0e-5,100


icalc = 1
sigma_c=0.1
sigma_d=0.1
for kk in range(len(Ra0_array)):
    Ra0 = Ra0_array[kk]
    f = f_array[kk]
    Deta = np.exp(Deta_array[kk])
    if normal:
        H = 0
    if internal:
        H = H_array[kk]
    cgeom = (1 + f + f**2)/3
    gamma=np.log(Deta)
    
    Ramed=Ra0*np.exp(0.5*gamma)
    #print(Ra0,Ramed,H,f,gamma,Rk)
    
    RT = Rk-1.0
    
    T_TDV = 1.0 - a*Rk**e/gamma/f**b
    
    
    Tm=T_TDV
    km = 1.0/(1.0 + RT*Tm)
    Ram = Ra0*np.exp(gamma*Tm)/km
    nite=0
    c=0
    d=0
    #print('---1---',Tm)
    if H>0:

       a1 = a1_Fpos
       a2 = a2_Fpos
       c = c_Fpos
       d = d_Fpos
       #print(a1,a2,c,d,icalc,itemax,nite,eps,f,cgeom,H,gamma,RT,Rk,Ra0,T_TDV,Tm,Ram,km)
       ite,Tm,Ram,test=solve_temperature(a1,a2,c,d,icalc,itemax,nite,eps,f,cgeom,H,gamma,RT,Rk,Ra0,T_TDV,Tm,Ram,km)
       #print(Tm,Ram,test)
       #print('---2---',Tm)
       sigma_a1 = sigma_a1_Fpos
       sigma_a2 = sigma_a2_Fpos
       sigma_c  = sigma_c_Fpos
       sigma_d  = sigma_d_Fpos

    aphi = aphi_Fpos
    bphi = bphi_Fpos
    cphi = cphi_Fpos
    dphi = dphi_Fpos
    ephi = ephi_Fpos
    Ftop = aphi*(Ram**bphi)*(f**dphi)/gamma**cphi/Rk**ephi
    Fbot = (Ftop - cgeom*H) / f**2.0
    
    sigma_aphi = sigma_aphi_Fpos
    sigma_bphi = sigma_bphi_Fpos
    sigma_cphi = sigma_cphi_Fpos
    sigma_dphi = sigma_dphi_Fpos
    sigma_ephi = sigma_ephi_Fpos
    
    if (H>0) and (Fbot<0): 
       nite = 0
       a1 = a1_Fneg
       a2 = a2_Fneg
       c = c_Fneg
       d = d_Fneg
       ite,Tm,Ram,test=solve_temperature(a1,a2,c,d,icalc,itemax,nite,eps,f,cgeom,H,gamma,RT,Rk,Ra0,T_TDV,Tm,Ram,km)           
       #print('---2---',ite,Tm,Ram,test)
       #print('---3---',Tm)
       aphi = aphi_Fneg
       bphi = bphi_Fneg
       cphi = cphi_Fneg
       dphi = dphi_Fneg
       ephi = ephi_Fneg
       Ftop = aphi*(Ram**bphi)*(f**dphi)/gamma**cphi/Rk**ephi
       Fbot = (Ftop - cgeom*H) / f**2

       sigma_aphi = sigma_aphi_Fneg
       sigma_bphi = sigma_bphi_Fneg
       sigma_cphi = sigma_cphi_Fneg
       sigma_dphi = sigma_dphi_Fneg
       sigma_ephi = sigma_ephi_Fneg

       sigma_a1 = sigma_a1_Fneg
       sigma_a2 = sigma_a2_Fneg
       sigma_c  = sigma_c_Fneg
       sigma_d  = sigma_d_Fneg

       if (Fbot>0):

           Fbot = 0.0
           Ftop = (1.0 + f + f**2)*H/3.0

           aphi = 0.5*(aphi_Fneg + aphi_Fpos)
           bphi = 0.5*(bphi_Fneg + bphi_Fpos)
           cphi = 0.5*(cphi_Fneg + cphi_Fpos)
           dphi = 0.5*(dphi_Fneg + dphi_Fpos)
           ephi = 0.5*(ephi_Fneg + ephi_Fpos)

           sigma_aphi = 0.5*(sigma_aphi_Fneg + sigma_aphi_Fpos)
           sigma_bphi = 0.5*(sigma_bphi_Fneg + sigma_bphi_Fpos)
           sigma_cphi = 0.5*(sigma_cphi_Fneg + sigma_cphi_Fpos)
           sigma_dphi = 0.5*(sigma_dphi_Fneg + sigma_dphi_Fpos)
           sigma_ephi = 0.5*(sigma_ephi_Fneg + sigma_ephi_Fpos)

           aux = Ftop*gamma**cphi/aphi/Ra0**bphi
           Tm  = np.log(aux)/gamma/bphi
    #print('---4---',Tm)
    #print(Ra0,',',H,',',f,',',format(Deta,'.1e'),',',1.0,Tm, format(Ram,'.3e'),km)
    err_a    = sigma_a/a
    err_c    = np.log(H)*sigma_c
    err_d    = np.log(Ram)*sigma_d
    err_aphi = sigma_aphi/aphi
    err_bphi = np.log(Ram)*sigma_bphi
    if (gamma > 0.0):
       err_cphi = abs(np.log(gamma))*sigma_cphi
    else:
       err_cphi = 0.0
    if (f > 0.0):
       err_b    = abs(np.log(f))*sigma_b
       err_dphi = abs(np.log(f))*sigma_dphi
    else:
       err_dphi = 0.0
    if (Rk > 0.0):
       err_e    = abs(np.log(Rk))*sigma_e
       err_ephi = abs(np.log(Rk))*sigma_ephi
    else:
       err_ephi = 0.0
    
    
    aux         = (cgeom*H)**c/Ram**d
    sigma_T_TDV = err_a + err_b + err_e

    #print( '----------> error: ',sigma_T_TDV)

    if (H > 0.0): 
       sigma_Tm = sigma_T_TDV + aux*(sigma_a1 + sigma_a2*f + (a1 + a2*f)*(err_c + err_d))
    else:
       sigma_Tm = sigma_T_TDV
    
    sigma_F = Ftop*(err_aphi + err_bphi + err_cphi + err_dphi + err_ephi)

    Urey = cgeom*H/Ftop

    if (Rk==1.0):
       Fcond_top = f + (f + 2.0)*H/6.0
    else:
       Fcond_top = np.log(1.0 + RT)/RT
    Nu = Ftop/Fcond_top
    
    print('---------->', kk, '<----------')
    print('Tm = ',round(Tm,5),' +/- ',round(sigma_Tm,5), Tm_array[kk])
    print('Ram = ',format(Ram,'.5e'))
    print('Ftop = ',round(Ftop,5),' +/- ',round(sigma_F,3),Ftop_array[kk])
    print('Fbot = ',round(Fbot,3) )
    print('Urey = ',round(Urey,5))
    print('                              ')
    
    ax.errorbar(Tm_array[kk], Tm, yerr=sigma_Tm, fmt='o', capsize=5, label='Data with Error Bars',color = 'k')
    ax2.errorbar(Ftop_array[kk],Ftop,yerr=sigma_F, fmt='o', capsize=5, label='Data with Error Bars',color = 'k')
    if (H>0):
        ax3.errorbar(Tm_array[kk], Tm, yerr=sigma_Tm, fmt='o', capsize=5,color = 'orange',ecolor='#849DAB',label='mix heat' if kk==12 else '')
    else:
        ax3.errorbar(Tm_array[kk], Tm, yerr=sigma_Tm, fmt='o', capsize=5,color = 'k',ecolor='#849DAB',label='basal heat' if kk==0 else '')


bwith=3
labelsize=25
for aa in [ax,ax2,ax3]:
    aa.tick_params(labelsize=labelsize)
    aa.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    # aa.grid()
ax.set_xlim(0.8,1.07)
ax.set_ylim(0.8,1.07)
ax2.set_xlim(0,10)
ax2.set_ylim(0,10)
ax.set_xlabel('observed interior temperature',fontsize = labelsize)
ax.set_ylabel('modeled interior temperature',fontsize = labelsize)
ax2.set_xlabel('F_obseverd/f',fontsize = labelsize)
ax2.set_ylabel('F_pred/f',fontsize = labelsize)
ax3.set_xlim(0.8,1)
ax3.set_ylim(0.8,1)
ax3.set_xlabel('observed interior temperature',fontsize = labelsize)
ax3.set_ylabel('modeled interior temperature',fontsize = labelsize)
ax3.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
ax3.legend(fontsize = labelsize)
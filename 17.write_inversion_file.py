#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 17:03:23 2023

@author: chingchen
"""
import os
import pandas as pd
import numpy as np
write_inversion_nlg_inverse_3v3p_solve_Tm = 0
write_scaling_output = 0
write_nlgi_1v4p_Tm_lid_input = 1
write_Tm_plot = 0
write_plot_Ram = 0

header_list = ['rsurf','f','Ea','H','Tm','ftop','fbot','Ur','raeff','dlid','Tlid']
if write_inversion_nlg_inverse_3v3p_solve_Tm:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    path = workpath+'nlg_inverse_3v3p_solve_Tm/'
    name = 'inverse_input_H.dat'
    name = '1124_basal_heat.dat'
    #name = '1124_mix_heat.dat'
    file = open(path+name, 'w')
    
    
    #model_information = pd.read_csv(workpath+'table_normal_jiching.dat',sep='\\s+')
    model_information = pd.read_csv(workpath+'model_information_all.csv',sep=',') 
    #model_information = pd.read_csv(workpath+'model_information_v2013_normal.csv',sep=',') 
    print('3','2','3', file=file)
    print('0.1d-2','0.1d-2', file=file)
    print('1.000','2.000','0.500','1.000','0.000','2.000', file=file)
    print('1.0d-06','100','1', file=file)
    for jj in range(len(model_information)):
        
        model = model_information.model[jj]
        rsurf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        gamma = round(np.log(Ea),3)
        if H==0 and Tm<1:
            
            print(model,rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03')
            print(rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03',file=file)
        
    file.close()

if write_scaling_output:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    path = workpath+'model_output_data_H_pure_basal.dat'
    file = open(path, 'w')
    
    model_information = pd.read_csv(workpath+'model_information_all.csv',sep=',')  
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        Ea = format(Ea,'.1e')
        raeff = format(raeff,'.3e')
        print(model,rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid)
        print(model,rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid,file=file)
        
    file.close()
if write_plot_Ram:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    model_information = pd.read_csv(workpath+'model_information_v2013_normal.csv',sep=',') 
    name = 'plot_Ram_v2013_result.dat'
    file = open(workpath+name, 'w')
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        gamma = round(np.log(Ea))
        print(rsruf,f,H,gamma,Tm,ftop,file=file)
    file.close()
    model_information = pd.read_csv(workpath+'model_information_all.csv',sep=',') 
    name = 'plot_Ram_v2013_result_all.dat'
    file = open(workpath+name, 'w')
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        gamma = round(np.log(Ea),3)
        print(rsruf,f,H,gamma,Tm,ftop,file=file)
    file.close()
    
if write_nlgi_1v4p_Tm_lid_input:
    workpath='/Users/chingchen/Desktop/StagYY_Works/'
    file = open(workpath+'nlg_inverse_1v4p_Tm_lid/nlgi_1v4p_Tm_lid_jiching.dat','w')
    print('4','1','4', file=file)
    print('0.5d-2','0.5d-3', file=file)
    print('2.000','2.001','0.000','2.000','0.250','0.001','1.000','0.001', file=file)
    print('1.0d-06','100','1', file=file)
    scaling_data = pd.read_csv(workpath+'model_output_data.dat',sep='\\s+',header=None,names=header_list)
    for jj in range(len(scaling_data)):
        if scaling_data.Tm[jj] <1:
            print(round(scaling_data.rsurf[jj],3),scaling_data.f[jj],round(np.log(scaling_data.Ea[jj]),3),
                  scaling_data.H[jj],scaling_data.Tm[jj],'1e-05',file=file)
    file.close()
    #----------------make input file----------------------
#     file = open(workpath+'nlg_inverse_1v4p_Tm_lid/'+'nlg_inverse_1v4p_Tm_lid.in','w')
#     print('nlgi_1v4p_Tm_lid_jiching.dat',file=file)
#     print('nlgi_1v4p_Tm_lid_jiching.out',file=file)
#     file.close()
#     cmd = '''
# cd %(workpath)snlg_inverse_1v4p_Tm_lid/
# ./nlg_inverse_1v4p_Tm_lid-MH-YY < nlg_inverse_1v4p_Tm_lid.in
# grep Parametre nlgi_1v4p_Tm_lid_jiching.out > parameters_Tm
# awk -F+ '{print$1}'  parameters_Tm | awk -F'=  ' '{print$2}' >parameters
# awk -F'- ' '{print$2}'  parameters_Tm > sigmas
# rm parameters_Tm
#     '''%locals()
#     print(cmd)
#     os.system(cmd)

if write_Tm_plot:
    workpath='/Users/chingchen/Desktop/StagYY_Works/'
    name = 'internal_for_Tm_jiching.dat'
    file = open(workpath+'nlg_inverse_3v3p_solve_Tm/'+name,'w')
    scaling_data = pd.read_csv(workpath+'model_output_data.dat',sep='\\s+',header=None,names=header_list)
    for jj in range(len(scaling_data)):
        print(round(scaling_data.rsruf[jj],3),round(np.log(scaling_data.Ea[jj]),3),scaling_data.f[jj],
              scaling_data.H[jj],scaling_data.Tm[jj],scaling_data.ftop[jj],scaling_data.fbot[jj],file=file)
    file.close()
    #Ra0_array, Deta_array, f_array, H_array,Tm_array, Ftop_array = np.loadtxt(path+file).T
    
    

# workpath = '/Users/chingchen/Desktop/StagYY_Works/'
# datapath = '/Users/chingchen/Desktop/data/'
# path = workpath+'nlg_inverse_3v3p_solve_Tm/'
# name = 'inverse_input_H.dat'
# name = '1118_basal_mix_heat.dat'
# file = open(path+name, 'w')



# header_list = ['model','rsurf','f','Ea','H','Tm','ftop','fbot','Ur','raeff','dlid','Tlid']
# model_information = pd.read_csv(path+'1118_model_data_H_pure_basal.dat',sep='\\s+',header=None,names=header_list)

# print('3','2','3', file=file)
# print('0.1d-2','0.1d-2', file=file)
# print('1.000','2.000','0.500','1.000','0.000','2.000', file=file)
# print('1.0d-06','100','1', file=file)
# for jj in range(len(model_information)):
    
#     model = model_information.model[jj]
#     Ea = model_information.Ea[jj]
#     rsurf = model_information.rsurf[jj]
#     f = model_information.f[jj]
#     Tm = model_information.Tm[jj]
#     ftop = model_information.ftop[jj]
#     fbot = model_information.fbot[jj]
#     H = model_information.H[jj]
#     gamma = round(np.log(Ea),3)
#     #rsurf = round(rsurf,3)
#     if Tm>0:
        
#         print(model,rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03')
#         #print(path+name)
#         print(rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03',file=file)
    
# file.close()
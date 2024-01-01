#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 17:03:23 2023

@author: chingchen
"""
# import os
import pandas as pd
import numpy as np
write_inversion_nlg_inverse_3v3p_solve_Tm = 0 # write input files for inversion
write_scaling_output = 0                      # write all output for models
write_plot_Ram = 0                            # 
write_nlgi_1v4p_Tm_lid_input = 0              # only for internal heat models
write_3v3p_solve_Tm_plot = 1                  # 

if write_inversion_nlg_inverse_3v3p_solve_Tm:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    datapath='/Users/chingchen/Desktop/StagYY_Works/data_scaling/'
    path = workpath+'nlg_inverse_3v3p_solve_Tm/'
    name = '1231_mix_heat.dat'
    file = open(path+name, 'w')
    model_information = pd.read_csv(workpath+'model_information_23end.csv',sep=',') 
    print('3','2','3', file=file)
    print('0.1d-2','0.1d-2', file=file)
    print('1.000','2.000','0.500','1.000','0.000','2.000', file=file)
    print('1.0d-06','100','1', file=file)
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsurf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid,tw1,tw2 = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        gamma = round(np.log(Ea),3)
        if  Tm<1:
            
            print(model,rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03')
            print(rsurf,gamma,f,Tm,'1e-05',ftop,'1e-03',file=file)
        
    file.close()

if write_scaling_output:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    datapath='/Users/chingchen/Desktop/StagYY_Works/data_scaling/'
    path = workpath+'1231_model_output_data_H_pure_basal.dat'
    file = open(path, 'w')
    
    model_information = pd.read_csv(workpath+'model_information_23end.csv',sep=',')  
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid,tw1,tw2 = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        Ea = format(Ea,'.1e')
        raeff = format(raeff,'.3e')
        #print(model,rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid)
        print(model,rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid,file=file)
        
    file.close()
if write_plot_Ram:
    workpath = '/Users/chingchen/Desktop/StagYY_Works/'
    datapath = '/Users/chingchen/Desktop/data/'
    datapath='/Users/chingchen/Desktop/StagYY_Works/data_scaling/'
    model_information = pd.read_csv(workpath+'model_information_23end.csv',sep=',') 
    name = 'plot_Ram.dat'
    file = open(workpath+name, 'w')
    for jj in range(len(model_information)):
        model = model_information.model[jj]
        rsruf,f,Ea,H,Tm,ftop,fbot,Ur,raeff,dlid,Tlid,tw1,tw2 = np.loadtxt(datapath+model+'_scaling_grep_data.txt')
        gamma = round(np.log(Ea),3)
        print(rsruf,f,H,gamma,Tm,ftop,file=file)
    file.close()
    
if write_nlgi_1v4p_Tm_lid_input: # only for internal heat models
    header_list = ['rsurf','f','Ea','H','Tm','ftop','fbot','Ur','raeff','dlid','Tlid']
    workpath='/Users/chingchen/Desktop/StagYY_Works/'
    name = '1231_mix_heat_1v4p_Tm.dat'
    file = open(workpath+'nlg_inverse_1v4p_Tm_lid/'+name,'w')
    print('4','1','4', file=file)
    print('0.5d-2','0.5d-3', file=file)
    print('2.000','2.001','0.000','2.000','0.250','0.001','1.000','0.001', file=file)
    print('1.0d-06','100','1', file=file)
    scaling_data = pd.read_csv(workpath+'1231_model_output_data_H_pure_basal.dat',
                               sep='\\s+',header=None,names=header_list)
    for jj in range(len(scaling_data)):
        if scaling_data.Tm[jj] <1 and scaling_data.H[jj] > 0:
            print(round(scaling_data.rsurf[jj],3),scaling_data.f[jj],
                  round(np.log(scaling_data.Ea[jj]),3),
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

if write_3v3p_solve_Tm_plot:
    header_list = ['rsurf','f','Ea','H','Tm','ftop','fbot','Ur','raeff','dlid','Tlid']
    workpath='/Users/chingchen/Desktop/StagYY_Works/'
    name = 'internal_for_Tm_jiching.dat'
    file = open(workpath+'nlg_inverse_3v3p_solve_Tm/'+name,'w')
    
    scaling_data = pd.read_csv(workpath+'1231_model_output_data_H_pure_basal.dat',
                               sep='\\s+',header=None,names=header_list)
    for jj in range(len(scaling_data)):
        print(round(scaling_data.rsurf[jj],3),round(np.log(scaling_data.Ea[jj]),3),
              scaling_data.f[jj],scaling_data.H[jj],scaling_data.Tm[jj],
              scaling_data.ftop[jj],scaling_data.fbot[jj],file=file)
        model = scaling_data.index[jj]
        name = 'solve_Tm_TDV-H_jiching_'+str(model)
        file2 = open(workpath+'nlg_inverse_3v3p_solve_Tm/'+name,'w')
        print('1', file=file2)
        print(round(scaling_data.rsurf[jj],3),',',scaling_data.H[jj],
              ',',scaling_data.f[jj],',',scaling_data.Ea[jj],',1.0',file=file2)
        print('1.23,1.5,0.1',file=file2)
        print('0.05,0.02,0.02',file=file2)
        print('6.0332,5.1235,1.00,0.250',file=file2)
        print('0.5141,0.6952,0.001,0.001',file=file2)
        print('4.36,3.00,1.72,0.333',file=file2)
        print('0.15,0.15,0.0,0.0',file=file2)
        print('1.934,0.339,1.87,0.0,0.82',file=file2) # parameters of ftop scaling 
        print('0.06,0.004,0.03,0.0,0.01',file=file2)
        print('1.57,0.270,1.21,0.0,0.82',file=file2)
        print('0.06,0.004,0.03,0.0,0.01',file=file2)
        print('1.0d-5,100',file=file2)
        file2.close()
    file.close()
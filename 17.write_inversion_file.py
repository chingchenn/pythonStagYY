#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 17:03:23 2023

@author: chingchen
"""

import pandas as pd
workpath = '/Users/chingchen/Desktop/StagYY_Works/nlgi_module/'

path = workpath+'normal_inversion.dat'
file = open(path, 'w')


model_output = pd.read_csv(workpath+'table_normal_jiching.dat',sep='\\s+')
print('3','2','3', file=file)
print('0.5d-3','0.5d-3', file=file)
print('1.000','2.000','0.500','1.000','0.000','2.000', file=file)
print('1.0d-06','100','1', file=file)
for jj in range(len(model_output)):
    
# start = time.time()
    f = model_output.f[jj]
    gamma = model_output.gamma[jj]
    rsurf = model_output.rsurf[jj]
    tmm = model_output.Tmm[jj]
    ftop = model_output.Ftop[jj]
    print(rsurf,gamma,f,tmm,'1e-05',ftop,'1e-03',file=file)
    
file.close()
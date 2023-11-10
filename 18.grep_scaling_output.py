#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 15:50:21 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
import sys
path = '/Users/chingchen/Desktop/data/'
path='/lfs/jiching/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/'
modelpath = '/Users/chingchen/Desktop/model/'
modelpath='/lfs/jiching/ScalingLaw_model/23agu/'
model= sys.argv[1]
file = path+model+'_time.dat'
ff=np.loadtxt(file,skiprows=1)
header_list=['istep','time','F_top','F_bot','Tmin',
             'Tmean','Tmax','Vmin','Vrms','Vmax',
             'eta_min','eta_mean',' eta_max','ra_eff',
             'Nu_top','Nu_bot','H_int','NON1','NON2',
             'NON3','NON4','NON5','NON6','NON7','NON8',
             'NON9','NON10']
data = np.loadtxt(file,skiprows=1)
kk = np.hsplit(data, 27)
mm=np.array(kk[0:16]).reshape(16,len(data))
ff1 = pd.DataFrame(mm.T, columns = header_list[0:16])
ff17 = pd.DataFrame(np.array(kk[-1]), columns = ['H_int'])
ff=pd.concat([ff1.istep,ff1.time,ff1.F_top, ff.F_bot,ff.Tmean,ff.Nu_top,ff.Nu_bot,ff1.ra_eff, ff17],axis=1)
ff.to_csv(path+model+'_scaling_time_data.csv')  

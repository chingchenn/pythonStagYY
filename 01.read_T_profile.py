#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:00:42 2023

@author: chingchen
"""

import pandas as pd
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

### PATH ###
local = 1
if local:
    path = '/Users/chingchen/Desktop/data/'
    workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution/'
    modelpath = '/Users/chingchen/Desktop/model/'
    figpath = '/Users/chingchen/Desktop/figure/'
else:
    path = '/lfs/jiching/data/'
    workpath = '/lfs/jiching/ScalingLaw_model/'
    modelpath = '/lfs/jiching/ScalingLaw_model/'
    figpath = '/lfs/jiching/figure/'

labelsize = 30
bwith = 3

fig_temperature_model = 0
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
              '#AE6378','#282130','#7E9680','#24788F',
              '#849DAB','#EA5E51','#35838D','#4198B9',
              '#414F67','#97795D','#6B0D47','#A80359','#52254F']
header_list = ['depth','T']
time_list = [0.000   0.225   0.453   0.680   0.908   1.135   1.363   1.590   1.818   2.045   2.273   2.500   2.728   2.955   3.183   3.410   3.638   3.865   4.093   4.320   4.548   4.550]
# read data file
data = pd.read_csv(workpath+'test_v1_T-profiles.dat',
    header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
data = data.replace('D','e',regex=True).astype(float) # data convert to float

fig,(ax,ax2) = plt.subplots(1,2,figsize=(12,12))


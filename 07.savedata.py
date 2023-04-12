#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:38:36 2023

@author: chingchen
"""

import stagpy
import os,sys
import numpy as np
from stagpy import stagyydata
from stagpy import field


model = 'w0210'
model = sys.argv[1]
path = '/lfs/jiching/'
savepath = '/lfs/jiching/data/'
data = stagyydata.StagyyData(path+model)
for shot in range(1,53):
    print(shot)

    kk1,kk2,kk3,kk4 = field.get_meshes_fld(data.snaps[shot],'T')
    eta1,eta2,eta3,eta4 = field.get_meshes_fld(data.snaps[shot],'eta')
    rho1,rho2,rho3,rho4 = field.get_meshes_fld(data.snaps[shot],'rho')
    np.savetxt(savepath+'/'+model+'_get_meshes_x_of_'+str(shot)+'.txt',kk1)
    np.savetxt(savepath+'/'+model+'_get_meshes_z_of_'+str(shot)+'.txt',kk2)
    np.savetxt(savepath+'/'+model+'_get_meshes_T_of_'+str(shot)+'.txt',kk3)
    np.savetxt(savepath+'/'+model+'_get_meshes_eta_of_'+str(shot)+'.txt',eta3)
    np.savetxt(savepath+'/'+model+'_get_meshes_rho_of_'+str(shot)+'.txt',rho3)

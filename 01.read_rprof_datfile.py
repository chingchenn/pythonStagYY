#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 22:44:59 2023

@author: chingchen
"""


import os
datapath = '/Users/chingchen/Desktop/data/'
path = '/Users/chingchen/Desktop/model/'
model = 'w0209'
for mm in range(1,2):
    model = 'w020'+str(mm)
    model = 'w0209'
    file = datapath+model+'_rprof.dat'
    sourceFileName = file
    sourceFileData = open(file,'r')
    ListOfLine = sourceFileData.read().splitlines()
    new_data = ListOfLine[1:]
    n = len(new_data)
    p = 129
    m = int(n/p)
    for i in range(m):
        destFileName = path+model+'/datafile/'
        if not os.path.isdir(path+model):
            os.mkdir(path+model)
        if not os.path.isdir(destFileName):
            os.mkdir(destFileName)
        destFileData = open(destFileName+model+'_data_'+str(i)+'.txt','w')
        
        if (i==m-1):
            for line in new_data[i*p+1:]:
                destFileData.write(line+'\n')
        else:
            for line in new_data[i*p+1:(i+1)*p]:
                destFileData.write(line+'\n')
        destFileData.close()
    print(model+'==DONE==')


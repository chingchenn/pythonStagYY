#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 22:44:59 2023

@author: chingchen
"""


import os, sys 
datapath = '/Users/chingchen/Desktop/data/'
#datapath = '/lfs/jiching/data/'
path = '/Users/chingchen/Desktop/model/'
#path = '/lfs/jiching/ScalingLaw_model/'

model = sys.argv[1]
p = sys.argv[2]
p = int(p) + 1 
#for mm in range(1,2):
    
    # model = 'w080'+str(mm)
    # model='ca06'
    # p = 129
    #model = 'h06_192'
    #p = 193
    # model = 'h04_256'
    # p = 257
file = datapath+model+'_rprof.dat'
sourceFileName = file
sourceFileData = open(file,'r')
ListOfLine = sourceFileData.read().splitlines()
new_data = ListOfLine[1:]
n = len(new_data)
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
print(model+'==DONE==',m-1)

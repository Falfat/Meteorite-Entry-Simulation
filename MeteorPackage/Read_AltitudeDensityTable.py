# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:30:18 2018

@author: Yujie
"""

import numpy as np

#function
def read(filename):
    file=open(filename,'r')
    z_list=[]
    pa_list=[]
    H_list=[]
    for line in file:
        words=line.split()
        try:
            if words[0]=='#':
                raise ValueError
            z=float(words[0])
            pa=float(words[1])
            H=float(words[2])
            z_list.append(z)
            pa_list.append(pa)
            H_list.append(H)
        except:
            pass
    file.close()
    return z_list,pa_list,H_list

def readcurvedata(filename):
    file=open(filename,'r')
    alt_list=[]
    energy_list=[]
    for line in file:
#        print("dingdong")
        words=line.split(",")
        try:
#            print(words[0])
            alt = float(words[0])
#            print(alt)
            energy = float(words[1])
            alt_list.append(alt)
            energy_list.append(energy)
        except:
            pass
    file.close()
    return np.array(alt_list), np.array(energy_list)

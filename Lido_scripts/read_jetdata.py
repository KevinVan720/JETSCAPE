# -*- coding: utf-8 -*-
"""
Created on Sat May 19 23:04:55 2018

@author: kevin
"""

import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):
    with open(filename) as f:
        rst=[]
        for line in f:
            strlist=line.split(sep=" ")
            rst.append([float(a) for a in strlist])
        rst=np.array(rst)
        return rst
            

data0=readfile("JShqjet_7000_0.5_0_0.5.dat")
data1=readfile("JShqjet_7000_0.5_0.5_1.dat")
data2=readfile("JShqjet_7000_0.5_1_1.5.dat")
data3=readfile("JShqjet_7000_0.5_1.5_2.dat")

fig, axs=plt.subplots(nrows=1, ncols=1, sharex=True)
axs.errorbar(data0[:,0], data0[:,1]*625, yerr=data0[:,2]*625, label="|y|<0.5 (*625)")
axs.errorbar(data1[:,0], data1[:,1]*125, yerr=data1[:,2]*125, label="0.5<|y|<1 (*125)")
axs.errorbar(data2[:,0], data2[:,1]*25, yerr=data2[:,2]*25, label="1<|y|<1.5 (*25)")
axs.errorbar(data3[:,0], data3[:,1]*5, yerr=data3[:,2]*5, label="1.5<|y|<2 (*5)")
axs.set_yscale('log')
axs.set_xscale('log')
axs.set_aspect(0.14)
axs.legend()
fig.savefig('cs.png')
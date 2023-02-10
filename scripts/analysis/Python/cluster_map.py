#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:56:54 2023

@author: jholt
"""
import numpy as np
import matplotlib.pylab as plt

Row={}
Row[1]=np.array([6,5,4,2,3]).astype(int)-1
Row[2]=np.array([14,23,15,22,1]).astype(int)-1
Row[3]=np.array([11,16,21,17,7]).astype(int)-1
Row[4]=np.array([12,20,18,19,8,9,10]).astype(int)-1
Row[5]=np.array([13]).astype(int)-1

plt.figure(figsize=[11.69,8.27])
gap=0.05
Position=np.zeros((23,4))

for ir in range(1,6):
    A=np.load('bounds.npz')
    bounds=A['arr_0']

    bnds=bounds[Row[ir],:]
    nP=np.size(Row[ir])
    if ir != 5:
        WW=np.sum(bnds[:,2])+(nP-1)*gap+2*gap
    gap_s=gap/WW

    x0=gap_s
    if ir ==1 :
        y0 = 1-x0*1.5

    Hmax=0.
    for ip in Row[ir]:
        W=bounds[ip,2]/WW
        H=W*bounds[ip,3]/bounds[ip,2]
        y1=y0-H
        Position[ip,0] = x0
        Position[ip,1] = y1
        Position[ip,2] = W
        Position[ip,3] = H    
        x0=x0+gap_s+W
        ax=plt.axes()
        ax.set_position(Position[ip,:])
        Hmax=np.max([H,Hmax])
    print(y0)    
    y0=y0 - Hmax -gap_s*1.5
np.savez('Position.npz',Position)
    
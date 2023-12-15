#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:14:27 2022

@author: willy
"""
import os
import numpy as np
direc = os.getcwd()
Asp = np.array([6.505,6.605,6.705,6.805,6.905,7.005,7.105,7.205,7.305,7.405,7.505,7.605,7.705,7.805,7.905,8.005,8.105,8.205,8.305,8.405,8.505])


base1=np.genfromtxt(open(direc+"/zcdom127.dfo", "r"))
base2=np.genfromtxt(open(direc+"/zcdom148.dfo", "r"))
base3=np.genfromtxt(open(direc+"/zcdom169.dfo", "r"))
base4=np.genfromtxt(open(direc+"/zcdom190.dfo", "r"))
base5=np.genfromtxt(open(direc+"/zcdom211.dfo", "r"))
base6=np.genfromtxt(open(direc+"/zcdom232.dfo", "r"))

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def encontrar(base1):
    Asp1 = []
    for j in range(len(Asp)):
        aux = find_nearest(base1[:,1],Asp[j])
        aux1 = np.where(base1[:,1]==aux)[0][0]
        print(base1[:,0][aux1])
        Asp1.append(base1[:,0][aux1])
        
print('casoRH1_C1')
encontrar(base1)
print('casoRH2_C1')
encontrar(base2)
print('casoRH3_C1')
encontrar(base3)
print('casoRH1_C2')
encontrar(base4)
print('casoRH2_C2')
encontrar(base5)
print('casoRH3_C2')
encontrar(base6)


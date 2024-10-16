# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 13:13:55 2024

@author: jackl
"""

#import sys
import numpy as np

def distance(x1,y1,z1,x2,y2,z2):
    dist=np.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2)
    return dist

with open('h_360.inp','r') as inp:
    f0=inp.readlines()

with open('formed_bchl361-79002.g96','r') as g96:
    g0=g96.readlines()

atoms=[]
for line in g0:
    if '360 BCL' in line:
        atoms.append(line)

def compare(inp):
    inp_x=float(inp.split()[2])/10
    inp_y=float(inp.split()[3])/10
    for line in atoms:
        x=float(line.split()[4])
        y=float(line.split()[5])
        if(x-inp_x > -0.0001 and x-inp_x < 0.0001) and (y-inp_y > -0.0001 and y-inp_y < 0.0001):
            #print(line.replace('\n',''))
            print(line.split()[2])

start=0
for line in f0:
    if '$end' in line:
        if(start==1):
            break
    elif(start==1):
        compare(line)
    elif 'C1' in line:
        start=1
        print('begin')
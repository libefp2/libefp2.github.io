# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:20:41 2024

@author: jackl
"""

import numpy as np
import sys

site='359'

def distance(x1,y1,z1,x2,y2,z2):
    dist=np.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2)
    return dist

def dist_check(residue):
    for atom in ring_atoms:
        for atom2 in residue:
            if(distance(atom2['x'],atom2['y'],atom2['z'],atom['x'],atom['y'],atom['z'])<1.5):
                return 1
    return 0

Rings=['MG','CHA','CHB','HB','CHC','HC','CHD','HD','NA','C1A','C2A','H2A','C3A','H3A','C4A',
       'CMA','HMA1','HMA2','HMA3','NB','C1B','C2B','C3B','C4B','CMB','HMB1','HMB2','HMB3',
       'CAB','OBB','CBB','HBB1','HBB2','HBB3','NC','C1C','C2C','H2C','C3C','H3C','C4C',
       'CMC','HMC1','HMC2','HMC3','CAC','HAC1','HAC2','CBC','HBC1','HBC2','HBC3','ND',
       'C1D','C2D','C3D','C4D','CMD','HMD1','HMD2','HMD3','CAD','OBD','CBD','HBD','CGD',
       'O1D','O2D','CED','HED1','HED2','HED3']

#with open('shell_bchl361-79002.gro','r') as inp:
#    f0=inp.readlines()

with open('bchl359-50028.g96','r') as g96:
    g0=g96.readlines()


#Get BCL atoms and Head ring atoms
bcl_atoms=[]
ring_atoms=[]
for line in g0:
    if(line.split()[0]==site):
        parts=line.split()
        atom_info = {
        'residue_number': int(parts[0]),
        'residue': parts[1],
        'atom_name': parts[2],
        'atom_number': int(parts[3]),
        'x': float(parts[4]),
        'y': float(parts[5]),
        'z': float(parts[6])}
        bcl_atoms.append(atom_info)
        if atom_info['atom_name'] in Rings:
            ring_atoms.append(atom_info)

#Get residues <15 Angstroms from the head ring
efp_res=[]
resi=[]
for line in g0:
    parts = line.split()
    if 'Generated' in line:
        continue
    elif(len(parts))>= 7:
        #print(line)
        if(int(parts[0])!=atom_info['residue_number']):
            if(dist_check(resi)):
                efp_res.append(resi[0])
                resi=[]
        atom_info = {
        'residue_number': int(parts[0]),
        'residue': parts[1],
        'atom_name': parts[2],
        'atom_number': int(parts[3]),
        'x': float(parts[4]),
        'y': float(parts[5]),
        'z': float(parts[6])}
        resi.append(atom_info)
        
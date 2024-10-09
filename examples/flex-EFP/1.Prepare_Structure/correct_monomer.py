# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:11:23 2023

@author: jackl
"""

#Script to add non-backbone components from OPM file to the CHARMGUI output
#CHARM has added lipid bilayer, but could not handle CLA, PQN, etc.
#OPM original file has the information for these, but they are translated in X and Y directions

#Need to move X: +12.802
#Need to move Y: -52.107

pdb='opm_5oy0.pdb'

f=open(pdb,'r')
f0=f.readlines()
f.close()

lipids='monomer_lipid.pdb'

g=open(pdb,'r')
g0=g.readlines()
g.close()

monomer='ABCDEFGHIJKLM'
start=0
addres=[]

for line in f0:
    if line[0:3]=='HET':
        start=1
        print(line)
        
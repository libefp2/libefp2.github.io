# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 07:42:51 2024

@author: jackl

Sample execution:
    cut_qm.py ala_33_473.inp a0001.efp

This file reads in an input file and a reference efp file that should be a good geometric match (low RMSD), 
AND the .inp file has listed atoms that should be deleted, or polarizable points that should be deleted.
Atoms that will be added to the QM region should be removed from the fragment file. Additionally, a bridge
between a QM-atom and an MM-atom that are covalently bound should have the QM atom removed, and some 
information of the MM-atom should be removed. The coordinates, monopole, and screen parameter sections 
for the MM-atom are retained, while the dipoles, quadrupoles, octupoles, and polarizable point sections 
should be removed.
"""

import sys
import numpy as np

def distance(x1,y1,z1,x2,y2,z2):
    dist=np.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2)
    return dist

def get_specials(lines): # Return the atom indexes of atoms to be removed entirely and atoms to remove polarizations
    atoms=[]
    pols=[]
    for line in lines:
        j=-1
        if 'erased:' in line:
            while(line.split()[j]!='erased:'):
                atoms.append(line.split()[j])
                j-=1
        elif 'remove:' in line:
            while(line.split()[j]!='remove:'):
                pols.append(line.split()[j])
                j-=1
    i=0
    atomid=[]
    polid=[]
    start=0
    for line in lines:
        if '$end' in line:
            start=0
        elif(start==1):
            i+=1
            if line.split()[0] in atoms:
                atomid.append(i)
            elif line.split()[0] in pols:
                polid.append(i)
            elif 'H000' in line:
                atomid.append(i)
        elif('C1'==line.split()[0]):
            start=1
    return atomid,polid

def get_coords(lines,atoms,pols):
    start=0
    coords=[]
    rem_coords=[]
    for line in lines:
        if 'STOP' in line:
            return coords,rem_coords
        elif(start==1) and (line[0]=='A'):
            atomnum=int(line[1:3])
            if atomnum in atoms:
                rem_coords.append(line)
            elif atomnum in pols:
                coords.append(line)
                rem_coords.append(line)
            else:
                coords.append(line)
        elif(start==1) and (line[0]=='B'):
            if(len(line.split()[0])==4):
                atomnum=int(line[2])
                atomnum2=int(line[3])
                if atomnum in atoms:
                    rem_coords.append(line)
                elif atomnum2 in atoms:
                    rem_coords.append(line)
                elif atomnum in pols:
                    rem_coords.append(line)
                    coords.append(line)
                elif atomnum2 in pols:
                    rem_coords.append(line)
                    coords.append(line)
                else:
                    coords.append(line)
            elif(len(line.split()[0])==5):
                atomnum=int(line[2:4])
                atomnum2=int(line[4])
                if atomnum in atoms:
                    rem_coords.append(line)
                elif atomnum2 in atoms:
                    rem_coords.append(line)
                elif atomnum in pols:
                    rem_coords.append(line)
                    coords.append(line)
                elif atomnum2 in pols:
                    rem_coords.append(line)
                    coords.append(line)
                else:
                    coords.append(line)
            elif(len(line.split()[0])==6):
                atomnum=int(line[2:4])
                atomnum2=int(line[4:6])
                if atomnum in atoms:
                    rem_coords.append(line)
                elif atomnum2 in atoms:
                    rem_coords.append(line)
                elif atomnum in pols:
                    rem_coords.append(line)
                    coords.append(line)
                elif atomnum2 in pols:
                    rem_coords.append(line)
                    coords.append(line)
                else:
                    coords.append(line)
        elif 'COORDINATES' in line:
            start=1

def get_monopoles(lines,coords):
    monopoles=[]
    keep_names=[]
    start=0
    for atom in coords:
        keep_names.append(atom.split()[0])
    for line in lines:
        if(start==1):
            if 'STOP' in line:
                return monopoles
            if line.split()[0] in keep_names:
                monopoles.append(line)
        if 'MONOPOLES' in line:
            start=1
            
def get_dipoles(lines,coords):
    dipoles=[]
    cut_names=[]
    start=0
    for atom in coords:
        cut_names.append(atom.split()[0])
    for line in lines:
        if(start==1):
            if 'STOP' in line:
                return dipoles
            if line.split()[0] in cut_names:
                continue
            else:
                dipoles.append(line)
        if 'DIPOLES' in line:
            start=1

def get_quadrupoles(lines,coords):
    quadrupoles=[]
    cut_names=[]
    start=0
    k=0
    j=0
    for atom in coords:
        cut_names.append(atom.split()[0])
    for line in lines:
        if(start==1):
            if 'STOP' in line:
                return quadrupoles
            elif(k==1):
                k-=1
            elif(j==1):
                quadrupoles.append(line)
                j-=1
            elif line.split()[0] in cut_names:
                k=1
            else:
                quadrupoles.append(line)
                j=1
        if 'QUADRUPOLES' in line:
            start=1

def get_octupoles(lines,coords):
    octupoles=[]
    cut_names=[]
    start=0
    j=0
    k=0
    for atom in coords:
        cut_names.append(atom.split()[0])
    for line in lines:
        if(start==1):
            if 'STOP' in line:
                return octupoles
            elif(k>0):
                k-=1
            elif(j>0):
                octupoles.append(line)
                j-=1
            elif line.split()[0] in cut_names:
                k=2
            else:
                octupoles.append(line)
                j=2
        if 'OCTUPOLES' in line:
            start=1

def get_polarpts(lines,cut_coords):
    polars=[]
    remove_names=[]
    start=0
    j=0
    for atom in cut_coords:
        remove_names.append(atom.split()[0])
    for line in lines:
        mindist=20.0
        if(start==1):
            if 'STOP' in line:
                return polars
            elif(j>0):
                polars.append(line)
                j-=1
            elif(line[0:2]=='CT'):
                for atom in cut_coords:
                    current_dist=distance(line.split()[1],line.split()[2],line.split()[3],atom.split()[1],atom.split()[2],atom.split()[3])
                    if(current_dist<mindist):
                        mindist=current_dist
                if(mindist>3):
                    polars.append(line)
                    j=3
        if 'POLARIZABLE POINTS' in line:
            start=1

def get_screen(lines,coords,title):
    params=[]
    keep_names=[]
    start=0
    for atom in coords:
        keep_names.append(atom.split()[0])
    for line in lines:
        if(start==1):
            if 'STOP' in line:
                return params
            if line.split()[0] in keep_names:
                params.append(line)
        if title in line:
            start=1

def get_header(lines):
    header=[]
    for line in lines:
        header.append(line)
        #print(line)
        if 'COORDINATES (BOHR)' in line:
            return header

def main(inp, efp):
    with open(inp,'r') as inp:
        inp_lines=inp.readlines()
    with open(efp,'r') as efp:
        efp_lines=efp.readlines()
    header=get_header(efp_lines)
    rem_atoms, rem_pols = get_specials(inp_lines)
    keep_coords,rem_coords=get_coords(efp_lines,rem_atoms,rem_pols)
    keep_monop=get_monopoles(efp_lines,keep_coords)
    keep_dip=get_dipoles(efp_lines,rem_coords)
    keep_quadrup=get_quadrupoles(efp_lines,rem_coords)
    keep_octup=get_octupoles(efp_lines,rem_coords)
    keep_pols=get_polarpts(efp_lines,rem_coords)
    keep_screen=get_screen(efp_lines,keep_coords,'SCREEN ')
    keep_screen2=get_screen(efp_lines,keep_coords,'SCREEN2')
    with open('cut_'+efp,'w') as outfile:
        for outline in header:
            outfile.write(outline)
        for outline in keep_coords:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      ' MONOPOLES\n')
        for outline in keep_monop:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      ' DIPOLES\n')
        for outline in keep_dip:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      ' QUADRPOLES\n')
        for outline in keep_quadrup:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      ' OCTUPOLES\n')
        for outline in keep_octup:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      ' POLARIZABLE POINTS\n')
        for outline in keep_pols:
            outfile.write(outline)
        outfile.write(' STOP\n'+
                      'SCREEN2      (FROM VDWSCL=   0.700)\n')
        for outline in keep_screen:
            outfile.write(outline)
        outfile.write('STOP\n'+
                      'SCREEN       (FROM VDWSCL=   0.700)\n')
        for outline in keep_screen2:
            outfile.write(outline)
        outfile.write('STOP\n'+
                      ' $END\n')
    
if __name__ == "__main__":
    #main()
    main(sys.argv[1], sys.argv[2])
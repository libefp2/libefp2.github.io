# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 2024

@author: jackl
"""
import sys
import numpy as np
import os

directory = '/depot/lslipche/data/yb_boss/flexible_efp/efpdb/'
ang_cutoff=0.20                     # Minimum RMSD allowed for "good" match
cutoff=ang_cutoff*1.8897259886      # Cutoff convert to Bohr

inp=sys.argv[1]    
#inp='g_952.inp'                 # a_22_304.inp for example. Input file with fragment coords
with open(inp, "r") as orig:
    orgs = orig.readlines()
    
outname=inp.replace('.inp','.efp')  # Output file name, same as input but changed extension.
                                    # If no match is found, then outname will not be used.

def get_RMSD(coords1,coords2,mass):
    tot=0.0
    length=len(coords1)
    dim=len(coords1[0])
    totmass=0.0
    for i in range(length):
        for j in range(dim):
            tot+=mass[i]*(coords1[i][j]-coords2[i][j])**2
        totmass+=mass[i]
    rmsd=np.sqrt(tot)/np.sqrt(totmass)
    return rmsd

def kabsch_algorithm(coords, coords2):
        
    com1 = np.mean(coords, axis=0)
    com2 = np.mean(coords2, axis=0)
    
    coords1_aligned = coords - com1
    coords2_aligned=coords2 - com2
    
    covariance_matrix = np.dot(coords1_aligned.T, coords2_aligned)
    #print(covariance_matrix)
    
    u, _, vh = np.linalg.svd(covariance_matrix)
    rotation_matrix = np.dot(u, vh)

    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix[:, -1] *= -1
    return rotation_matrix, com1, com2

def apply_transform(coords, rotation_matrix, com1, com2):
    translated_coords = coords - com1

    rotated_coords = np.dot(translated_coords, rotation_matrix)
    rotated_coords+=com2
    return rotated_coords

amino_acid_dict = {'a':'ala','r':'arg','n':'asn','d':'asp','c':'cys',
                   'q':'gln','e':'glu','g':'gly','h':'hip','i':'ile',
                   'l':'leu','k':'lys','m':'met','f':'phe','p':'pro',
                   's':'ser','t':'thr','w':'trp','y':'tyr','v':'val',
                   'hp':'hip','hd':'hid','he':'hie'}
atom_weights = {'H':1.0,'O':15.9949100,'C':12.0000000,'N':12.0030700,'S':32.0650000,'M':23.3040000}

try:
    res=amino_acid_dict[inp.split('_')[0]]
except:
    print('No library folder for: '+inp)
    exit()

directory+=res

file_extension = '.efp'  # change this to your desired file extension


#Grab coordinates from the input file, convert from Angstroms to bohrs.
start=0
xyz=[]
orig_coords=[]
weights=[]
inp_at_lines=[]
for line in orgs:
    if '$end' in line:
        start=0
    elif(start==1):
        inp_at_lines.append(line)
        xyz.append(float(line.split()[2])/0.529177249)
        xyz.append(float(line.split()[3])/0.529177249)
        xyz.append(float(line.split()[4])/0.529177249)
        orig_coords.append(xyz)
        xyz=[]
        weights.append(atom_weights[line.split()[0][0]])
    elif 'C1' in line:
        start=1
minrmsd=20.0     # Arbitrary starting value. This will be replaced with real RMSD values in the loop.

# Loop through files in the directory
for filename in os.listdir(directory):
    if filename.endswith(file_extension):
        full_path = os.path.join(directory, filename)
        with open(full_path, "r") as term:
            terms = term.readlines()
    else:
        continue
    start=0
    xyz=[]
    term_coords=[]
    for line in terms:
        if 'BO21' in line:
            start=0
        elif(start==1):
            xyz.append(float(line.split()[1]))
            xyz.append(float(line.split()[2]))
            xyz.append(float(line.split()[3]))
            term_coords.append(xyz)
            xyz=[]
        elif 'COORDINATES (BOHR)' in line:
            start=1
    new_coords=np.array(term_coords)
    coords2=np.array(orig_coords)
    rotation_matrix, com1, com2 = kabsch_algorithm(new_coords, coords2)
    aligned_new_coords = apply_transform(new_coords, rotation_matrix, com1,com2)
    rmsd=(get_RMSD(aligned_new_coords,coords2,weights))

    if(rmsd<minrmsd):
        minrmsd=rmsd 
        match=full_path

if(minrmsd<cutoff):
    os.system("python step4.Flexible_V5.py "+full_path+' '+inp+' -d -xr')
else:
    print('no match, run GAMESS for: '+inp)
    #os.system('gms_slurm -p 20 -q standby -v 2023 '+inp)

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

inp=sys.argv[1]                     # a_22_304.inp for example. Input file with fragment coords
with open(inp, "r") as orig:
    orgs = orig.readlines()
    
outname=inp.replace('.inp','.efp')  # Output file name, same as input but changed extension.
                                    # If no match is found, then outname will not be used.

def rotate_vector(vector, rotation_matrix):
    #Different from coords because dipoles assumed centered on coordinate. NOT placed in Cartesian space
    return np.dot(vector,rotation_matrix)

def rotate_quadrupole(quadrupole, rotation_matrix):
    #R x Q x R_t
    #Quadrupoles are 3x3s. But .efp give only six values, order is assumed xx, xy, xz, yy, yz, zz
    #With symmetric off-diagonal elements.
    q = np.array([
        [quadrupole[0], quadrupole[3], quadrupole[4]],
        [quadrupole[3], quadrupole[1], quadrupole[5]],
        [quadrupole[4], quadrupole[5], quadrupole[2]],
    ])
    q_rotated = rotation_matrix.T @ q @ rotation_matrix
    return [q_rotated[0, 0], q_rotated[1, 1], q_rotated[2, 2],
            q_rotated[0, 1], q_rotated[0, 2], q_rotated[1, 2]]

def rotate_octupole(octupole, rotation_matrix):
    #Oxyz=Rxi*Ryj*Rzk*Oijk summed over i=x,y,z;j=x,y,z;z=x,y,z;
    #Each term Oxyz of new tensor is composed of 
    #27 pieces generated from old tensor and rotation matrix!
    #The indexes are hard-coded for some better readiability
    o = np.array([
        [[octupole[0], octupole[3], octupole[4]], [octupole[3], octupole[5], octupole[9]], [octupole[4], octupole[9], octupole[7]]],
        [[octupole[3], octupole[5], octupole[9]], [octupole[5], octupole[1], octupole[6]], [octupole[9], octupole[6], octupole[8]]],
        [[octupole[4], octupole[9], octupole[7]], [octupole[9], octupole[6], octupole[8]], [octupole[7], octupole[8], octupole[2]]],
    ])
    #fancy function below does what we want. Look up np.einsum if you need clarification.
    #o_rotated = np.einsum('ij,jkl,lm->ikm', rotation_matrix, o, rotation_matrix.T)
    o_rotated = np.zeros((3, 3, 3))
    for p in range(0,3):
        for q in range(0,3):
            for r in range(0,3):
                o_rotated[p][q][r]=single_term_oct(p,q,r,o,rotation_matrix.T)
    return [
        o_rotated[0, 0, 0], o_rotated[1, 1, 1], o_rotated[2, 2, 2],
        o_rotated[0, 0, 1], o_rotated[0, 0, 2], o_rotated[0, 1, 1],
        o_rotated[1, 1, 2], o_rotated[0, 2, 2], o_rotated[1, 2, 2],
        o_rotated[0, 1, 2]
    ]

def single_term_oct(p,q,r,octupole,rot_mat):
    term=0.0
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                term+=rot_mat[p][i]*rot_mat[q][j]*rot_mat[r][k]*octupole[i][j][k]
    return term

def rotate_polarizability(polarizability, rotation_matrix):
        #RxQxR_t
    q = np.array([
        [polarizability[0], polarizability[3], polarizability[4]],
        [polarizability[6], polarizability[1], polarizability[5]],
        [polarizability[7], polarizability[8], polarizability[2]]])
    q_rotated = rotation_matrix.T @ q @ rotation_matrix
    return [q_rotated[0, 0],q_rotated[1, 1],q_rotated[2, 2],
            q_rotated[0, 1],q_rotated[0, 2],q_rotated[1, 2],
            q_rotated[1, 0],q_rotated[2, 0],q_rotated[2, 1]]
    #return q_rotated.flatten()

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

    #if np.linalg.det(rotation_matrix) < 0:
    #    rotation_matrix[:, -1] *= -1
    return rotation_matrix, com1, com2

def apply_transform(coords, rotation_matrix, com1, com2):
    translated_coords = coords - com1

    rotated_coords = np.dot(translated_coords, rotation_matrix)
    rotated_coords+=com2
    return rotated_coords

amino_acid_dict = {'a':'ala','r':'arg','n':'asn','d':'asp','c':'cys',
                   'q':'gln','e':'glu','g':'gly','h':'his','i':'ile',
                   'l':'leu','k':'lys','m':'met','f':'phe','p':'pro',
                   's':'ser','t':'thr','w':'trp','y':'tyr','v':'val',
                   'hp':'hip','hd':'hid','he':'hie'}
atom_weights = {'H':1.0,'O':15.9949100,'C':12.0000000,'N':12.0030700,'S':32.0650000,'M':23.3040000}

try:
    res=amino_acid_dict[inp[0]]
except:
    print('No library folder for: '+inp)
#directory = '/depot/lslipche/data/yb_boss/flexible_efp/efpdb/'+res
directory+=res
file_extension = '.efp'  # change this to your desired file extension

start=0
xyz=[]
orig_coords=[]
weights=[]
for line in orgs:
    if '$end' in line:
        start=0
    elif(start==1):
        xyz.append(float(line.split()[2])/0.529177249)
        xyz.append(float(line.split()[3])/0.529177249)
        xyz.append(float(line.split()[4])/0.529177249)
        orig_coords.append(xyz)
        xyz=[]
        weights.append(atom_weights[line.split()[0][0]])
    elif 'C1' in line:
        start=1

minrmsd=20.0

# Loop through files in the directory
for filename in os.listdir(directory):
    if filename.endswith(file_extension):
        full_path = os.path.join(directory, filename)
        #full_path=directory+'/gly5817.efp'
        #print(f"Processing file: {full_path}")
        with open(full_path, "r") as term:
            terms = term.readlines()

    start=0
    xyz=[]
    term_coords=[]
    for line in terms:
        if 'BO21' in line:
            start=0
        elif(start==1):
            #print(line.split()[1])
            #print(float(line.split()[1])*0.529177249)
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
    #print(rmsd)
    if(rmsd<minrmsd):
        minrmsd=rmsd
        if(rmsd<cutoff):    
            match=full_path
            #break
    #print(minrmsd)
#%%

#If good match found, transform that .efp to coords in .inp
# Coords, dipoles, quadrupoles, octupoles, polarizable points and tensor
# Monopoles and screen/screen2 do not need changed

#print(rmsd/1.8897259886, full_path, inp)
#print(get_RMSD(aligned_new_coords,coords2,weights))
#print(aligned_new_coords)
#print(coords2)

if(rmsd<cutoff):
    with open(match, "r") as efp:
        f0 = efp.readlines()
    start=0
    #coords=[]
    temp=[]
    for line in f0:
        newline=line
        if ' STOP' in line:
            start=0
        elif(start==1):
            #coords.append(apply_transform([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],rotation_matrix,com1,com2))
            coords=(apply_transform([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],rotation_matrix,com1,com2))
            newline=line.split()[0]
            while(len(newline)<7):
                newline=newline+' '
            for i in range(len(coords)):
                newline+=str(f"{coords[i]:16.10f}")
                #newline=newline.replace(line.split()[i+1],str(f"{coords[i]:16.10f}"))
            newline+=str(f"{float(line.split()[4]):12.7f}")
            newline+=str(f"{float(line.split()[5]):5.1f}")+'\n'
        elif(start==2):
            dipole=rotate_vector([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],rotation_matrix)
            newline=line.split()[0]
            while(len(newline)<7):
                newline=newline+' '
            for i in range(len(dipole)):
                newline+=f"{dipole[i]:16.10f}"
            newline+='\n'
        elif(start==3):
            newline=''
            if '>' in line:
                name=line.split()[0]
                while(len(name)<7):
                    name=name+' '
                quadrupole=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])]
            else:
                quadrupole.append(float(line.split()[0]))
                quadrupole.append(float(line.split()[1]))
                rotated_quad=rotate_quadrupole(quadrupole, rotation_matrix)
                temp.append(name+'   '+f"{rotated_quad[0]:14.10f}"+'  '+f"{rotated_quad[1]:14.10f}"+'  '
                            +f"{rotated_quad[2]:14.10f}"+'  '+f"{rotated_quad[3]:14.10f}"+' >\n'+
                            '          '+f"{rotated_quad[4]:14.10f}"+'  '+f"{rotated_quad[5]:14.10f}"+'\n')
        elif(start==4):
            newline=''
            if(j%3==0):
                name=line.split()[0]
                while(len(name)<7):
                    name=name+' '
                #print(line)
                octupole=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])]
                j+=1
            elif(j%3==1):
                for i in range(4):
                    octupole.append(float(line.split()[i]))
                j+=1
            else:
                octupole.append(float(line.split()[0]))
                octupole.append(float(line.split()[1]))
                rotated_oct=rotate_octupole(octupole, rotation_matrix)
                #print(rotated_oct)
                temp.append(name+'     '+f"{rotated_oct[0]:13.9f}"+'    '+f"{rotated_oct[1]:13.9f}"+'    '
                            +f"{rotated_oct[2]:13.9f}"+'    '+f"{rotated_oct[3]:13.9f}"+' >\n'+
                            '            '+f"{rotated_oct[4]:13.9f}"+'    '+f"{rotated_oct[5]:13.9f}"'    '
                            +f"{rotated_oct[6]:13.9f}"+'    '+f"{rotated_oct[7]:13.9f}"+' >\n'+
                            '            '+f"{rotated_oct[8]:13.9f}"+'    '+f"{rotated_oct[9]:13.9f}"+'\n')
                j=0
        elif(start==5):
            if 'CT' in line:
                newline=line.split()[0]
                while(len(newline)<4):
                    newline=newline+' '
                coords=(apply_transform([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],rotation_matrix,com1,com2))
                for i in range(len(coords)):
                    newline+=f"{coords[i]:16.10f}"
                newline+='\n'
                temp.append(newline)
            elif(j%3==0):
                polar=[float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
                j+=1
            elif(j%3==1):
                for i in range(4):
                    polar.append(float(line.split()[i]))
                j+=1
            else:
                polar.append(float(line.split()[0]))
                rotated_polar=rotate_polarizability(polar,rotation_matrix)
                temp.append('  '+f"{rotated_polar[0]:14.10f}"+'  '+f"{rotated_polar[1]:14.10f}"+'  '+
                            f"{rotated_polar[2]:14.10f}"+'  '+f"{rotated_polar[3]:14.10f}"+' >\n'+
                            '  '+f"{rotated_polar[4]:14.10f}"+'  '+f"{rotated_polar[5]:14.10f}"+'  '+
                            f"{rotated_polar[6]:14.10f}"+'  '+f"{rotated_polar[7]:14.10f}"+' >\n'+
                            f"{rotated_polar[8]:16.10f}"+'\n')
                j=0
        elif 'COORDINATES (BOHR)' in line:
            start=1
        elif ' DIPOLES' in line:
            start=2
        elif ' QUADRUPOLES' in line:
            start=3
            temp.append(line)
        elif ' OCTUPOLES' in line:
            temp.append(line)
            start=4
            j=0
        elif ' POLARIZABLE POINTS' in line:
            start=5
            temp.append(line)
            j=0
        if(start<3):
            temp.append(newline)
    
    outfile=inp.split('.')[0]+'.efp'
    with open(outfile,'w') as outfile:
        for line1 in temp:
            outfile.write(line1)
else:
    with open('no_matches.txt','a') as outfile2:
        outfile2.write(inp+'\n')

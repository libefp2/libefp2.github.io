# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:13:19 2024

@author: jackl

Creates .inp files of amino acids or molecules in the EFP region. Non-amino acids
molecules are treated with no virtual bonds (no broken bonds).

"""

import numpy as np
import sys

#efp_g96='efp_pair53004.g96'
efp_g96=sys.argv[1]
with open(efp_g96, 'r') as inp:
    f0 = inp.readlines()

#full_g96='confout_pair53004.g96'
full_g96=sys.argv[2]
with open(full_g96, 'r') as g96:
    g0 = g96.readlines()

#settings_file='2user_defined.txt'
settings_file=sys.argv[3]

#topol_file='edit_topol.itp'
topol_file=sys.argv[4]

# Data and parameters

#Charges of residues. NVAL and CGLN should be redundant, now. New logic assumes and N-terminal is +1 and
# any C-terminal is -1. spec_AAs are residues that have non-zero charge, AA-charge is the dictionary with
#that non-zero charge.
AA_charge = {'ASP': '-1', 'GLU': '-1', 'NVAL': '1', 'HIP': '1', 'LYS': '1', 'ARG': '1', 'CGLN': '-1'}
spec_AAs = ['ASP', 'GLU', 'NVAL', 'HIP', 'LYS', 'ARG', 'CGLN']

amino_acid_dict = {
    'ALA': 'a', 'ARG': 'r', 'ASN': 'n', 'ASP': 'd', 'CYS': 'c', 'GLN': 'q', 'GLU': 'e', 'GLY': 'g', 'HIS': 'h', 
    'ILE': 'i', 'LEU': 'l', 'LYS': 'k', 'MET': 'm', 'PHE': 'f', 'PRO': 'p', 'SER': 's', 'THR': 't', 'TRP': 'w', 
    'TYR': 'y', 'VAL': 'v', 'HIP': 'hp', 'HID': 'hd', 'HIE': 'he', 'NVAL': 'v', 'CGLN': 'q'
}

known_amino_acids = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'HIP', 'HID', 'HIE'
]

at_sym = {'H': '1.0', 'C': '6.0', 'N': '7.0', 'O': '8.0', 'MG': '12.0', 'P': '15.0', 'S': '16.0', 'FE': '26.0', 'NA': '11.0', 'CL': '17.0'}

def qm_atoms(set_file):
    qm=0
    QMs=[]
    pol_rem=[]
    border=0
    i=3
    with open(set_file,'r') as input_file:
        inp_lines=input_file.readlines()
    for line in inp_lines:
        if(border==1):
            if 'boundary' in line:
                i=0
            elif(i<1):
                if int(line.split()[3]) in QMs:
                    i+=1
                else:
                    print('Atom number: '+line.split()[3]+' not found in QM atoms. Check bridge atom order.')
                    break
            elif(i==1):
                QMs.append(int(line.split()[3]))
                pol_rem.append(int(line.split()[3]))
        elif(qm==1):
            if 'QM-MM' in line:
                border=1
                qm=0
            if(len(line.split())>3):
                QMs.append(int(line.split()[3]))
        elif 'QM_atoms' in line:
            qm=1
    return QMs, pol_rem

def cut_frag(head, tail):
    """Calculate virtual coordinates for a fragment based on desired bond distance."""
    desired_dist = 1.07886
    xh, yh, zh = [float(head.split()[i]) * 10 for i in range(4, 7)]
    xt, yt, zt = [float(tail.split()[i]) * 10 for i in range(4, 7)]
    dist_mag = np.sqrt((xh - xt)**2 + (yh - yt)**2 + (zh - zt)**2)
    
    h_t = [
        ((xt - xh) * desired_dist / dist_mag) + xh,
        ((yt - yh) * desired_dist / dist_mag) + yh,
        ((zt - zh) * desired_dist / dist_mag) + zh
    ]
    return h_t

def make_inp(fragment, QMs, POLs):
    """Generate input file for a fragment."""
    txt = []
    num_virtuals=0
    for atom in fragment:
        if 'H000' in atom:
            num_virtuals+=1
    if(len(QMs)+num_virtuals==len(fragment)):
        return
    if fragment[4].split()[1] in spec_AAs:
        charge = AA_charge[fragment[4].split()[1]]
    elif len(fragment[4].split()[1])==4:
        if(fragment[4].split()[1][0]=='N'):
            charge=1
        elif(fragment[4].split()[1][0]=='C'):
            charge=-1
    else:
        charge = 0
    
    if fragment[4].split()[1] in known_amino_acids:
        filename = amino_acid_dict[fragment[4].split()[1]] + '_' + fragment[4].split()[0] + '_' + fragment[0].split()[3] + '.inp'
    else:
        filename = fragment[4].split()[1].lower() + '_' + fragment[4].split()[0] + '_' + fragment[0].split()[3] + '.inp'
    
    txt.append(f" $contrl units=angs local=boys runtyp=makefp \n"
               f"       mult=1 icharg={charge} coord=cart icut=11 $end\n"
               f" $system timlim=99999 mwords=200 $end\n"
               f" $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06 $end\n"
               f" $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n"
               f" $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n"
               f" $MAKEFP POL=.t. DISP=.f. CHTR=.f. EXREP=.f. $end\n"
               f" $data\n {filename.split('.')[0]}\n C1\n")
    
    for atom in fragment:
        if 'H000' in atom:
            txt.append(atom)
        else:
            col1 = f" {atom.split()[2]}"
            col1 = col1.ljust(6)
            
            if atom.split()[2][0] == 'M':
                col2 = at_sym['MG']
            else:
                col2 = at_sym[atom.split()[2][0]]
            
            x, y, z = [float(atom.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(17)
            col4 = f"{y:.8f}".rjust(18)
            col5 = f"{z:.8f}".rjust(18)
            
            txt.append(f"{col1}{col2}{col3}{col4}{col5}\n")
    
    txt.append(" $end \n!comment atoms to be erased:")
    for atomname in QMs:
        txt.append(" "+atomname)
    txt.append("\n")
    txt.append("!polarization points to remove:")
    for atomname in POLs:
        txt.append(" "+atomname)
        
               #" $end\n")
    
    with open(filename, 'w') as outfile:
        outfile.writelines(txt)

def QM_MM_covalent(MM,QMs,topol_file):
    with open(topol_file, 'r') as lines:
        bonds=0
        bond_IDs=[]
        for line in lines:
            if '[ bonds ]' in line:
                bonds=1
                continue
            elif bonds and len(line.split())<2:
                return bond_IDs
            elif bonds and (line[0]==';'):
                continue
            elif bonds and (int(line.split()[0])==MM) and int(line.split()[1]) not in QMs:
                bond_IDs.append(int(line.split()[1]))
            elif bonds and (int(line.split()[1])==MM) and int(line.split()[0]) not in QMs:
                bond_IDs.append(int(line.split()[0]))
                
            

#Gather from user which atoms are QM region and which are MM covalently bound to QM
qm_IDs, MM_remove = qm_atoms(settings_file)

#Take MM atoms that are covalently bound to QM region, then find the other MM atoms bound to that
#border atom (MM-MM-QM bonds, looking for the first MM atom). These atoms will have polarization points 
#removed. The middle MM atom will be removed entirely to avoid overlap with QM virtual hydrogens.
pol_remove=[]
for atom in MM_remove:
    found_pol_remove=QM_MM_covalent(atom,qm_IDs,topol_file)
    for bound_atom in found_pol_remove:
        pol_remove.append(bound_atom)

# Process EFP shell residue IDs and names
efp_resis = []
prevres = '0'
start = 0

for line in f0:
    if start == 1:
        if(len(line.split())<3):
            break
        if(line.split()[0]!=prevres):
            efp_resis.append(line.split()[0])
            prevres = line.split()[0]
    if 'POSITION' in line:
        start = 1

# Fragment processing and input generation
i = 0
prev_co = []
frag = []
CAs = []
end = 0
start = 0

for line in g0:
    if 'END' in line:
        if(start==1):
            break
    elif 'SOL' in line or 'QSL' in line:
        continue
    elif start == 1:
        if line.split()[2] == 'CA':
            CAs.append(line)
        elif line.split()[2] == 'O':
            prev_co.append(line)
        elif line.split()[2] == 'C':
            prev_co.append(line)
            if(line.split()[0]==efp_resis[i]):
                if len(line.split()[1]) == 4:
                    if line.split()[1][1:4] in known_amino_acids:
                        if line.split()[1][0] == 'N':  # N-terminal residue
                            vH2 = cut_frag(CAs[-1], line)
                            frag.append(f" H000    1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
                        elif line.split()[1][0] == 'C':  # C-terminal residue
                            vH1 = cut_frag(frag[0], CAs[-2])
                            frag.append(f" H000    1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
                elif line.split()[1] in known_amino_acids:
                    #print(frag)
                    vH1 = cut_frag(frag[0], CAs[-2])
                    vH2 = cut_frag(CAs[-1], line)
                    frag.append(f" H000 1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
                    frag.append(f" H000 1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
                    
                if len(frag) > 1:
                    qm_names=[]
                    pol_names=[]
                    for atom in frag:
                        if 'H000' in atom:
                            continue
                        elif int(atom.split()[3]) in qm_IDs:
                            qm_names.append(atom.split()[2])
                        elif int(atom.split()[3]) in pol_remove:
                            pol_names.append(atom.split()[2])
                    make_inp(frag,qm_names,pol_names)
                    frag = []
                    i += 1
        if(len(frag)>1) and (frag[-1].split()[0])==efp_resis[i]:
            if(frag[0].split()[0]!=line.split()[0]) and frag[0].split()[1] not in known_amino_acids:
                if frag[0].split()[1][1:] not in known_amino_acids:
                    qm_names=[]
                    pol_names=[]
                    for atom in frag:
                        if int(atom.split()[3]) in qm_IDs:
                            qm_names.append(atom.split()[2])
                        elif int(atom.split()[3]) in pol_remove:
                            pol_names.append(atom.split()[2])
                    make_inp(frag,qm_names,pol_names)
                    frag = []
                    i += 1

        if line.split()[0] == efp_resis[i]:
            if len(prev_co) > 1 and line.split()[1] in known_amino_acids:
                #frag.extend(prev_co)
                frag.append(prev_co[-2])
                frag.append(prev_co[-1])
                prev_co = []
            frag.append(line)
    elif 'POSITION' in line:
        start = 1

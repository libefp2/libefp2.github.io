# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 12:28:23 2024

@author: jackl
"""

import sys

water_and_ions={'SOL':-0.834,'CL':-1.0,'NA':1.0}

def get_EFPs(efplines):
    efp_resis=[]
    start=0
    prevres='0'
    for line in efplines:
        if 'END' in line:
            if(len(efp_resis)>1):
                return efp_resis
        elif start == 1:
            if(line.split()[0]!=prevres):
                resid=line.split()[0]
                resname=line.split()[1]
                efp_resis.append([resid,resname])
                prevres = resid
        if 'POSITION' in line:
            start = 1

def get_MM_coords(mm_indexes,g96):
    MMs=[]
    temp_MMs=[]
    prev_MM=0
    i=0
    prev_resid=0
    for line in g96:
        if(len(line.split())<4) or (line[0]!=' '):
            continue
        elif(line.split()[0]!=prev_resid) and (prev_MM==1):
            i+=1
        if(line.split()[2]=='C') or (line.split()[2]=='O'):
            col1=line.split()[2][0]+line.split()[3]
            while(len(col1)<7):
                col1+=' '
            col2=float(line.split()[4])*18.897161646321
            col3=float(line.split()[5])*18.897161646321
            col4=float(line.split()[6])*18.897161646321
            temp_MMs.append(col1+'%20.12f'%col2+'%20.12f'%col3+'%20.12f'%col4+'\n')
        else:
            if(line.split()[0]==mm_indexes[i][0]) and (line.split()[1]==mm_indexes[i][1]):
                if(len(temp_MMs)>1):
                    MMs.append(temp_MMs[-2])
                    MMs.append(temp_MMs[-1])
                    temp_MMs=[]
                col1=line.split()[2][0]+line.split()[3]
                while(len(col1)<7):
                    col1+=' '
                col2=float(line.split()[4])*18.897161646321
                col3=float(line.split()[5])*18.897161646321
                col4=float(line.split()[6])*18.897161646321
                MMs.append(col1+'%20.12f'%col2+'%20.12f'%col3+'%20.12f'%col4+'\n')
                prev_MM=1
            else:
                prev_MM=0
        prev_resid=line.split()[0]
    return MMs


def get_MM_charges(mm_indexes,topol):
    MMs=[]
    temp_MMs=[]
    prev_MM=0
    i=0
    prev_resid=0
    for line in topol:
        if '[ bonds ]' in line:
            atomid=int(MMs[-1].split()[0][1:-1]+MMs[-1].split()[0][-1])
            while(i<len(mm_indexes)):
                if(mm_indexes[i][1]=='SOL'):
                    atomid+=1
                    col1=('O'+str(atomid))
                    while(len(col1)<7):
                        col1+=' '
                    col2=water_and_ions['SOL']
                    MMs.append(col1+'%21.10f\n' % col2)
                    col2=col2/(-2)
                    atomid+=1
                    col1=('H'+str(atomid))
                    while(len(col1)<7):
                        col1+=' '
                    MMs.append(col1+'%21.10f\n' % col2)
                    atomid+=1
                    col1=('H'+str(atomid))
                    while(len(col1)<7):
                        col1+=' '
                    MMs.append(col1+'%21.10f\n' % col2)
                else:
                    atomid+=1
                    col1=(mm_indexes[i][1][0]+str(atomid))
                    while(len(col1)<7):
                        col1+=' '
                    col2=water_and_ions[mm_indexes[i][1]]
                    MMs.append(col1+'%21.10f\n' % col2)
                i+=1
            return MMs
        if(line[0]!=' '):
            continue
        elif(line.split()[2]!=prev_resid) and (prev_MM==1):
            i+=1
            prev_MM=0
        if(line.split()[1]=='C') or (line.split()[1]=='O'):
            col1=line.split()[1][0]+line.split()[0]
            while(len(col1)<7):
                col1+=' '
            col2=float(line.split()[6])
            temp_MMs.append(col1+'%21.10f\n' % col2)
        else:
            if(line.split()[2]==mm_indexes[i][0]) and (line.split()[3]==mm_indexes[i][1]):
                if(len(temp_MMs)>1):
                    MMs.append(temp_MMs[-2])
                    MMs.append(temp_MMs[-1])
                    temp_MMs=[]
                col1=(line.split()[4][0]+line.split()[0])
                while(len(col1)<7):
                    col1+=' '
                col2=float(line.split()[6])
                MMs.append(col1+'%21.10f\n' % col2)
                prev_MM=1
            else:
                prev_MM=0
        prev_resid=line.split()[2]
    return MMs
               
def get_screen(charges):
    screens=[]
    for atom in charges:
        col1=atom.split()[0]
        while(len(col1)<7):
            col1+=' '
        col2='   1.0000000000  10.0000000000\n'
        screens.append(col1+col2)
    return screens
    
def main(efp_g96, full_g96, topol):
    with open(efp_g96,'r') as inp:
        shell_lines=inp.readlines()
    with open(full_g96,'r') as efp:
        full_lines=efp.readlines()
    with open(topol,'r') as top:
        topol_lines=top.readlines()
    
    efp_residues=get_EFPs(shell_lines)
    all_residues=get_EFPs(full_lines)
    mm_residues=[]
    for res in all_residues:
        if(res[1]=='XXX'):
            continue
        elif res in efp_residues:
            continue
        else:
            mm_residues.append(res)
    MM_coords=get_MM_coords(mm_residues,full_lines)
    MM_charge=get_MM_charges(mm_residues,topol_lines)
    MM_screen2=get_screen(MM_charge)
    with open('test_mm.efp','w') as outfile:
        outfile.write(' $PROT\n')
        outfile.write('TITLE\n')
        outfile.write('  COORDINATES (BOHR)\n')
        for outline in MM_coords:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write('MONOPOLES\n')
        for outline in MM_charge:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write('SCREEN2      (FROM VDWSCL=   0.700)\n')
        for outline in MM_screen2:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write(' $END')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    #main('efp_pair53004.g96','confout_pair53004.g96','edit_topol.itp')
    #    efp_structure, full_strucuter, topol.top
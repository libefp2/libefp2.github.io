# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:41:27 2024

@author: jackl
"""

import numpy as np
import sys

def cut_frag(head,tail):
    desired_dist=1.07886
    #print(head)
    xh=float(head.split()[4])*10
    yh=float(head.split()[5])*10
    zh=float(head.split()[6])*10
    xt=float(tail.split()[4])*10
    yt=float(tail.split()[5])*10
    zt=float(tail.split()[6])*10
    dist_mag=np.sqrt(((xh-xt)**2)+((yh-yt)**2)+((zh-zt)**2))
    #print(dist_mag,desired_dist)
    #t_h=[((xh-xt)*desired_dist/dist_mag)+xt,((yh-yt)*desired_dist/dist_mag)+yt,((zh-zt)*desired_dist/dist_mag)+zt]
    h_t=[((xt-xh)*desired_dist/dist_mag)+xh,((yt-yh)*desired_dist/dist_mag)+yh,((zt-zh)*desired_dist/dist_mag)+zh]
    return h_t

def make_inp(fragment):
    txt=[]
    if line.split()[0] in spec_AAs:
        charge=AA_charge(fragment[0].split()[1])
    else:
        charge='0'
    #print(line.split()[1])
    #print(amino_acid_dict[fragment[0].split()[1]])
    #filename=amino_acid_dict[fragment[0].split()[1]]+'_'+fragment[0].split()[0]+'.inp'
    filename=amino_acid_dict[line.split()[1]]+'_'+line.split()[0]+'.inp'
    txt.append(' $contrl units=angs local=boys runtyp=makefp \n'+
'       mult=1 icharg='+charge+' coord=cart icut=11 $end\n'+
' $system timlim=99999   mwords=200 $end\n'+
' $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06  $end\n'+
' $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n'+
' $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n'+
' $MAKEFP  POL=.t. DISP=.f. CHTR=.f.  EXREP=.f. $end\n'+
' $data\n'+
' '+filename.split('.')[0]+'\n'+
' C1\n')
    for atom in fragment:
        if 'H000' in atom:
            txt.append(atom)
        else:
            if(atom.split()[2][0]=='M'):
                col1=' '+atom.split()[2][0:2]+atom.split()[3]+' '
                col2=at_sym['MG']
                filename='h_'+atom.split()[0]+'.inp'
            else:
                col1=' '+atom.split()[2][0]+atom.split()[3]+'   '
                col2=at_sym[atom.split()[2][0]]
            x=float(atom.split()[4])*10
            #col3=(str(x))
            col3=str("{:.8f}".format(x))
            while(len(col3)<17):
                col3=' '+col3
            y=float(atom.split()[5])*10
            col4=str("{:.8f}".format(y))
            while(len(col4)<18):
                col4=' '+col4
            z=float(atom.split()[6])*10
            col5=str("{:.8f}".format(z))
            while(len(col5)<18):
                col5=' '+col5
            txt.append(col1+col2+col3+col4+col5+'\n')
    txt.append(' $end \n'+
' $comment Atoms to be erased:  $end\n')
    with open(filename,'w') as outfile:
        for line1 in txt:
            outfile.write(line1)

Rings=['MG','CHA','CHB','HB','CHC','HC','CHD','HD','NA','C1A','C2A','H2A','C3A','H3A','C4A',
       'CMA','HMA1','HMA2','HMA3','NB','C1B','C2B','C3B','C4B','CMB','HMB1','HMB2','HMB3',
       'CAB','OBB','CBB','HBB1','HBB2','HBB3','NC','C1C','C2C','H2C','C3C','H3C','C4C',
       'CMC','HMC1','HMC2','HMC3','CAC','HAC1','HAC2','CBC','HBC1','HBC2','HBC3','ND',
       'C1D','C2D','C3D','C4D','CMD','HMD1','HMD2','HMD3','CAD','OBD','CBD','HBD','CGD',
       'O1D','O2D','CED','HED1','HED2','HED3']

AA_charge={'ASP':'-1','GLU':'-1','NVAL':'1','HIP':'1','LYS':'1','ARG':'1','CGLN':'-1'}
spec_AAs=['ASP','GLU','NVAL','HIP','LYS','ARG','CGLN']
BCL_resids=['359','360','361','362','363','364','365','366',
            '725','726','727','728','729','730','731','732',
            '1091','1092','1093','1094','1095','1096','1097','1098']

amino_acid_dict = {'ALA':'a','ARG':'r','ASN':'n','ASP':'d','CYS':'c',
                   'GLN':'q','GLU':'e','GLY':'g','HIS':'h','ILE':'i',
                   'LEU':'l','LYS':'k','MET':'m','PHE':'f','PRO':'p',
                   'SER':'s','THR':'t','TRP':'w','TYR':'y','VAL':'v',
                   'HIP':'h','HID':'h','HIE':'h','NVAL':'v','CGLN':'q'}

at_sym={'H':'1.0','C':'6.0','N':'7.0','O':'8.0','MG':'12.0','P':'15.0','S':'16.0','FE':'26.0',
        'NA':'11.0','CL':'17.0'}

with open('editconf_bchl_361-79002.g96','r') as inp:
    f0=inp.readlines()

with open('bchl361-79002.g96','r') as g96:
    g0=g96.readlines()


#Get EFP shell resid and resnames
efp_resis=[]
prevres='0'
start=0
for line in f0:
    #print(line)
    if '13.32071' in line:
        break
    elif(start==1):
        #print(line)
        if(line.split()[0]!=prevres):
            resid=''
            resname=''
            for char in line.split()[0]:
                if(char.isalpha()):
                    resname+=char
                else:
                    resid+=char
                '''
                resi={
                    'residue_number': resid,
                    'residue': resname,}
                '''
            efp_resis.append(resid)
        prevres=line.split()[0]
    if 'POSITION' in line:
        start=1

i=0
prev_co=[]
frag=[]
prev_co=[]
CAs=[]
end=0
start=0
residues=[]
for line in g0:
    if 'END' in line:
        end+=1
        if(end==2):
            break
    #print(line)
#    if(line.split()[2]=='CA') or (line.split()[2]=='N'):
#
    if 'BCL' in line:
        continue
    elif 'SOL' in line:
        continue
    elif(start==1):
        #print(line)
        if(line.split()[2]=='CA'):
            CAs.append(line)
        elif(line.split()[2]=='O'):
            prev_co.append(line)
        #elif(line.split()[2]=='C') or (line.split()[2]=='O'):
        elif(line.split()[2]=='C'):
            prev_co.append(line)
            if line.split()[0] in efp_resis:
                print(line.split()[0])
                if(line.split()[1]=='NVAL'):
                    vH2=cut_frag(CAs[-1],line)
                    frag.append(' H000    1.0      '+str("{:.8f}".format(vH2[0]))+'       '+str("{:.8f}".format(vH2[1]))+'       '+str("{:.8f}".format(vH2[2]))+'\n')
                elif(line.split()[1]=='CGLN'):
                    vH1=cut_frag(frag[0],CAs[-2])
                    frag.append(' H000    1.0      '+str("{:.8f}".format(vH1[0]))+'       '+str("{:.8f}".format(vH1[1]))+'       '+str("{:.8f}".format(vH1[2]))+'\n')
                else:
                    #print(line)
                    vH1=cut_frag(frag[0],CAs[-2])
                    vH2=cut_frag(CAs[-1],line)
                    frag.append(' H000    1.0      '+str("{:.8f}".format(vH1[0]))+'       '+str("{:.8f}".format(vH1[1]))+'       '+str("{:.8f}".format(vH1[2]))+'\n')
                    frag.append(' H000    1.0      '+str("{:.8f}".format(vH2[0]))+'       '+str("{:.8f}".format(vH2[1]))+'       '+str("{:.8f}".format(vH2[2]))+'\n')


            if line.split()[0] in efp_resis:
            #if(line.split()[0]==efp_resis[i+1]):
                if(len(frag)>1):
                    make_inp(frag)
                    #residues.append(frag)
                    frag=[]
                    i+=1


        #elif(line.split()[0]==efp_resis[i]['residue_number']):
        #print(line.split()[0],efp_resis[i])
        #print(efp_resis[i])
        if efp_resis[i] in BCL_resids:
            i+=1
        if(line.split()[0]==efp_resis[i]):
            if(len(prev_co)>1):
                '''
                if 'NVAL' in line:
                    continue
                else:
                '''
                if(line.split()[1]!='NVAL'):
                    for atom in prev_co:
                        frag.append(atom)
                    prev_co=[]
            frag.append(line)
    elif 'POSITION' in line:
        start=1
        
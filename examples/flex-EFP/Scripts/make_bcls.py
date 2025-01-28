# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:42:41 2024

@author: jackl

reads in .g96 of efp region and creates .inp files for bacteriochlorophyl a molecules.
Head and tail groups are separate fragments.
"""

import numpy as np
import sys

#site='359'
g96_file=sys.argv[1]
#g96_file='efp_pair53004.g96'
#site=sys.argv[2]

tailside='C6'
headside='C5'
site=['752','754']
oldres='752'

def cut_frag(head,tail):
    desired_dist=1.07886
    xh=float(head.split()[4])*10
    yh=float(head.split()[5])*10
    zh=float(head.split()[6])*10
    xt=float(tail.split()[4])*10
    yt=float(tail.split()[5])*10
    zt=float(tail.split()[6])*10
    dist_mag=np.sqrt(((xh-xt)**2)+((yh-yt)**2)+((zh-zt)**2))
    #print(dist_mag,desired_dist)
    t_h=[((xh-xt)*desired_dist/dist_mag)+xt,((yh-yt)*desired_dist/dist_mag)+yt,((zh-zt)*desired_dist/dist_mag)+zt]
    h_t=[((xt-xh)*desired_dist/dist_mag)+xh,((yt-yh)*desired_dist/dist_mag)+yh,((zt-zh)*desired_dist/dist_mag)+zh]
    return h_t,t_h

def make_inp(fragment):
    txt=[]
    filename='cla_t_'+fragment[0].split()[0]+'.inp'
    txt.append(' $contrl units=angs local=boys runtyp=makefp \n'+
'       mult=1 icharg=0 coord=cart icut=11 $end\n'+
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
                filename='cla_h_'+atom.split()[0]+'.inp'
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
            
at_sym={'H':'1.0','C':'6.0','N':'7.0','O':'8.0','MG':'12.0','P':'15.0','S':'16.0','FE':'26.0',
        'NA':'11.0','CL':'17.0'
      }

Rings=['MG','CHA','CHB','CHC','CHD','NA','C1A','C2A','C3A','C4A','CMA','NB','C1B','C2B','C3B',
       'C4B','CMB','CAB','CBB','NC','C1C','C2C','C3C','C4C','CMC','CAC','CBC','ND','C1D','C2D',
       'C3D','C4D','CMD','CAD','OBD','CBD','CGD','O1D','O2D','CED','CAA','CBA','CGA','O1A',
       'O2A','C1','C2','C3','C4','C5','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11',
       'H12','H13','H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24','H25',
       'H26','H27','H28','H29','H30','H31','H32','H33','H34','H35','H36','H37','H38','H39',
       'H40','H41']

#with open('bchl359-50028.g96','r') as g96:
#    g0=g96.readlines()

with open(g96_file,'r') as g96:
    g0=g96.readlines()
    
#tailside and headside represent the atoms that are bound but are to be broken to create 
#separate fragment files.
curr_head=[]
curr_tail=[]
#BCLs=[]
for line in g0:
    if 'CLA' in line:
        if(oldres!=line.split()[0]):
            head_coord,tail_coord = cut_frag(head_cut,tail_cut)
            curr_head.append(' H000    1.0      '+str("{:.8f}".format(head_coord[0]))+'       '+str("{:.8f}".format(head_coord[1]))+'       '+str("{:.8f}".format(head_coord[2]))+'\n')
            curr_tail.append(' H000    1.0      '+str("{:.8f}".format(tail_coord[0]))+'       '+str("{:.8f}".format(tail_coord[1]))+'       '+str("{:.8f}".format(tail_coord[2]))+'\n')
            #if(oldres!=site):
            if oldres not in site:
                make_inp(curr_head)
                make_inp(curr_tail)
            curr_head=[]
            curr_tail=[]
            oldres=line.split()[0]
        atomname=line.split()[2]
        if atomname in Rings:
            curr_head.append(line)
        else:
            curr_tail.append(line)
        if(atomname==tailside):
            tail_cut=line
        elif(atomname==headside):
            head_cut=line
        
        
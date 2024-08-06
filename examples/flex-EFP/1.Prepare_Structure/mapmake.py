import numpy as np
import sys

g96=sys.argv[1]
site=sys.argv[2]

def dist(x1,x2,y1,y2,z1,z2):
    num = np.sqrt(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))
    return num

f = open(g96,'r')
f0 = f.readlines()

xbcl=[]
ybcl=[]
zbcl=[]
final=[]
start=0
i=0
oldnum=0
resume=0
names=['MG','CHA','CHB','HB','CHC','HC','CHD','HD','NA','C1A','C2A','H2A','C3A','H3A','C4A','CMA','HMA1','HMA2','HMA3','NB','C1B','C2B','C3B','C4B','CMB','HMB1','HMB2','HMB3','CAB','OBB','CBB','HBB1','HBB2','HBB3','NC','C1C','C2C','H2C','C3C','H3C','C4C','CMC','HMC1','HMC2','HMC3','CAC','HAC1','HAC2','CBC','HBC1','HBC2','HBC3','ND','C1D','C2D','C3D','C4D','CMD','HMD1','HMD2','HMD3','CAD','OBD','CBD','HBD','CGD','O1D','O2D','CED','HED1','HED2','HED3']


for line in f0:
    if(line.split()[0]==site):
        for atom in names:
            if(line.split()[2]==atom):
                xbcl.append(float(line.split()[4]))
                ybcl.append(float(line.split()[5]))
                zbcl.append(float(line.split()[6]))

for line in f0:
    if(oldnum!=line.split()[0]):
        resume=1
        oldnum=line.split()[0]
    if(line.split()[0]=='1099'):
        start=0
    if(line.split()[0]=='1'):
        start=1
    if(start==1 and line.split()[0]!=site and resume==1):
        while(i<len(xbcl)):
            if(dist(xbcl[i],float(line.split()[4]),ybcl[i],float(line.split()[5]),zbcl[i],float(line.split()[6]))<1.5):
                resname=line.split()[1][0]+'_'+line.split()[0]+'\n'
                if(line.split()[1]=='NVAL'):
                    resname='v_'+line.split()[0]+'\n'
                if(line.split()[1]=='CGLN'):
                    resname='g_'+line.split()[0]+'\n'
                if(line.split()[1]=='BCL'):
                    resname='h_'+line.split()[0]+'\n'+'name = t_'+line.split()[0]+'\n'
                resname='name = '+resname
                resname=resname.lower()
                final.append(resname)
                i=500
                resume=0
            i=i+1
    i=0

text=open('newlist.map','w')
for line in final:
    text.write(line)
text.close()
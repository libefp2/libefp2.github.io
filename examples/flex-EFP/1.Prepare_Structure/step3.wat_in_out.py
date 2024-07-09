import sys,math
'''
python step3.wat_in_out.py new_bchl359-71156.dat
output --> new_new_bchl359-71156.dat
'''

def count_ef_aas(pn):
	if pn == 359:
		return 155 
	elif pn == 360:
		return 179
	elif pn == 361:
		return 178
	elif pn == 362:
		return 196
	elif pn == 363:
		return 181
	elif pn == 364:
		return 176
	elif pn == 365:
		return 214
	elif pn == 366:
		return 138
	elif pn == 1098:
		return 138
	else:
		print('not available site', str(pn))
		exit()

def compute_distance(qmx,qmy,qmz,solv,pn):
	nfrag           = count_ef_aas(pn)
	efrag,dist,back = [],[],[]
	for i in range(0,len(solv),3):
		line  = solv[i]
		x,y,z = float(line.split()[4]),float(line.split()[5]),float(line.split()[6])
		temp  = []
		for j in range(len(qmx)):
			dx,dy,dz = x-qmx[j],y-qmy[j],z-qmz[j]
			r2       = pow(dx,2)+pow(dy,2)+pow(dz,2)
			temp.append(r2)
		if min(temp) <= 225:
			efrag.append(line.split()[3])
			dist.append(min(temp))
			back.append(min(temp))
		else:
			continue
	real_min = []
	for i in range(len(efrag)):
		if nfrag+i+1 >= 500:
			break
		else:
			n_wat = i+1
			real_min.append(efrag[back.index(min(back))])
			back[back.index(min(back))] = 1000000.0

	re_write = []
	for i in range(len(solv)):
		line = solv[i]
		if line.split()[3] in real_min:
			if line.split()[2] == 'QSL':
				if pn == 360:
					resn = 'QSL'
				else:
					resn = 'INS'
			else:
				resn = 'INS'
			x1 = '{:{align}{width}}'.format('%s'%line.split()[0],align='<',width=2)
			x2 = '{:{align}{width}}'.format('%s'%line.split()[1],align='>',width=5)
			x3 = '{:{align}{width}}'.format('%s'%resn,align='>',width=5)
			x4 = '{:{align}{width}}'.format('%s'%line.split()[3],align='>',width=7)
			x5 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[4])),align='>',width=15)
			x6 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[5])),align='>',width=15)
			x7 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[6])),align='>',width=15)
			x8 = '{:{align}{width}}'.format('%.10f'%float(line.split()[7]),align='>',width=15)
			re_write.append(x1+x2+x3+x4+x5+x6+x7+x8)
		else:
			re_write.append(line)
	return re_write,len(real_min),len(efrag),nfrag

def main(inputfile):
	ifile    = open(inputfile,'r')
	reflines = ifile.readlines()
	pn       = int(inputfile.split('-')[0].split('bchl')[1])

	if pn == 366:
		pn = 1098
	else:
		pn = pn

	bact_chlorin = ['MG','CHA','CHB','CHC','CHD','NA','C1A','C2A','C3A','C4A','CMA','NB','C1B','C2B','C3B','C4B','CMB','NC','C1C','C2C','C3C','C4C','CMC','CAC','CBC','ND','C1D','C2D','C3D','C4D','CMD','CAD','OBD','CBD','CGD','O1D','O2D','CED']

	qmx,qmy,qmz    = [],[],[]
	prot,solv,ions = [],[],[]
	for obj in reflines:
		line = obj.strip()
		if line.split()[2] == 'BCL':
			prot.append(line)
			if int(line.split()[3]) == pn:
				if line.split()[1] in bact_chlorin:
					qmx.append(float(line.split()[4]))
					qmy.append(float(line.split()[5]))
					qmz.append(float(line.split()[6]))
				else:
					continue
			else:
				continue
		elif line.split()[2] == 'QSL':
			solv.append(line)
		elif line.split()[2] == 'SOL':
			solv.append(line)
		elif line.split()[2] == 'CL':
			ions.append(line)
		else:
			prot.append(line)

	wats,n_wat,r_wat,n_aas = compute_distance(qmx,qmy,qmz,solv,pn)
	pdb                    = prot+wats+ions
	print(inputfile,'total efrag: '+str(n_aas)+' '+str(n_wat)+' '+str(n_aas+n_wat))
	print(inputfile,'cut waters : '+str(r_wat-n_wat)+'\n'+'-'*50)
	ofile = open('new_'+inputfile,'w')
	for i in range(len(pdb)):
		ofile.write(pdb[i]+'\n')

main(sys.argv[1])

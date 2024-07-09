import sys,math
'''
python step4.prepare_qmefp.py list359.map new_new_bchl359-71156.dat
output --> new_new_new_bchl359-71156.dat
'''

def wrt_env(line,env):
	x = '{:{align}{width}}'.format('%s'%env,align='>',width=5)
	return line+x

def rewrt_env(line,env):
	x1 = '{:{align}{width}}'.format('%s'%line.split()[0],align='<',width=2)
	x2 = '{:{align}{width}}'.format('%s'%line.split()[1],align='>',width=5)
	x3 = '{:{align}{width}}'.format('%s'%line.split()[2],align='>',width=5)
	x4 = '{:{align}{width}}'.format('%s'%line.split()[3],align='>',width=7)
	x5 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[4])),align='>',width=15)
	x6 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[5])),align='>',width=15)
	x7 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[6])),align='>',width=15)
	x8 = '{:{align}{width}}'.format('%.10f'%float(line.split()[7]),align='>',width=15)
	x9 = '{:{align}{width}}'.format('%s'%env,align='>',width=5)
	return x1+x2+x3+x4+x5+x6+x7+x8+x9

def chg_sat(line,chg):
	x1 = '{:{align}{width}}'.format('%s'%line.split()[0],align='<',width=2)
	x2 = '{:{align}{width}}'.format('%s'%line.split()[1],align='>',width=5)
	x3 = '{:{align}{width}}'.format('%s'%line.split()[2],align='>',width=5)
	x4 = '{:{align}{width}}'.format('%s'%line.split()[3],align='>',width=7)
	x5 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[4])),align='>',width=15)
	x6 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[5])),align='>',width=15)
	x7 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[6])),align='>',width=15)
	x8 = '{:{align}{width}}'.format('%.10f'%float(chg),align='>',width=15)
	x9 = '{:{align}{width}}'.format('%s'%line.split()[8],align='>',width=5)
	return x1+x2+x3+x4+x5+x6+x7+x8+x9

def main(inputfile1,inputfile2):
	ifile1    = open(inputfile1,'r') # efp list ==> .map
	ifile2    = open(inputfile2,'r') # dat
	reflines1 = ifile1.readlines()
	reflines2 = ifile2.readlines()
	pn        = inputfile2.split('-')[0].split('bchl')[1]
	if pn == '366':
		pn = '1098'
	else:
		pn = pn

	efrag = []
	for obj in reflines1:
		line = obj.strip().split()
		efrag.append(line[2])

	pdb  = []
	temp = []
	for obj in reflines2:
		line = obj.strip()
		if line.split()[2] == 'NVAL':
			res1 = 'v_'+line.split()[3]
			res2 = res1
		elif line.split()[2] == 'CGLN':
			res1 = 'g_'+line.split()[3]
			res2 = res1
		elif line.split()[2] == 'BCL':
			res1 = 'h_'+line.split()[3]
			res2 = 't_'+line.split()[3]
		elif line.split()[2] == 'INS' or line.split()[2] == 'SOL' or line.split()[2] == 'QSL' or line.split()[2] == 'CL':
			res1 = ''
			res2 = ''
		else:
			res1 = line.split()[2][0:1].lower()+'_'+line.split()[3]
			res2 = res1

		if line.split()[3] == pn:
			pdb.append(wrt_env(line,'QM'))
			if line.split()[0] == 'MG':
				x1,y1,z1 = float(line.split()[4]),float(line.split()[5]),float(line.split()[6])
		elif res1 in efrag or res2 in efrag:
			pdb.append(wrt_env(line,'EFP'))
		elif line.split()[2] == 'INS':
			pdb.append(wrt_env(line,'EFP'))
		else:
			if pn == '360' and line.split()[2] == 'QSL' and line.split()[3] == '1099':
				pdb.append(wrt_env(line,'QM'))
			#	if line.split()[0] == 'O':
			#		x2,y2,z2 = float(line.split()[4]),float(line.split()[5]),float(line.split()[6])
			#		dx,dy,dz = x2-x1,y2-y1,z2-z1
			#		r2       = dx**2+dy**2+dz**2
			#		if math.sqrt(r2) <= 2.5:
			#			pdb.append(wrt_env(line,'QM'))
			#			temp.append('qm')
			#		else:
			#			pdb.append(wrt_env(line,'EFP'))
			#			temp.append('efp')
			#	else:
			#		print inputfile2,temp[0]
			#		if 'qm' in temp:
			#			pdb.append(wrt_env(line,'QM'))
			#		else:
			#			pdb.append(wrt_env(line,'EFP'))
			else:
				pdb.append(wrt_env(line,'MM'))

#	print inputfile2+' System charge: '+str(round(tot2+mmchg,5))
	print(inputfile2)
	ofile = open('new_'+inputfile2,'w')
	for i in range(len(pdb)):
		ofile.write(pdb[i]+'\n')


main(sys.argv[1],sys.argv[2])

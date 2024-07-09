import sys
'''
example execution command
python step1.define_ff_chg.py topol.top bchl359-71156.g96 
output --> bchl359-71156.dat
'''

def data_format(line,chg):
	if line.split()[2] == 'MG':
		atom = 'MG'
	elif line.split()[2] == 'CL':
		atom = 'CL'
	else:
		atom = line.split()[2][0:1]

	x1 = '{:{align}{width}}'.format('%s'%atom,align='<',width=2)   
	x2 = '{:{align}{width}}'.format('%s'%line.split()[2],align='>',width=5)   
	x3 = '{:{align}{width}}'.format('%s'%line.split()[1],align='>',width=5)   
	x4 = '{:{align}{width}}'.format('%s'%line.split()[0],align='>',width=7)   
	x5 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[4])*10),align='>',width=15)
	x6 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[5])*10),align='>',width=15)
	x7 = '{:{align}{width}}'.format('%.8f'%(float(line.split()[6])*10),align='>',width=15)
	x8 = '{:{align}{width}}'.format('%.10f'%float(chg),align='>',width=15)

	return x1+x2+x3+x4+x5+x6+x7+x8+'\n'

def main(inputfile1,inputfile2):
	ifile1    = open(inputfile1,'r') # topology
	ifile2    = open(inputfile2,'r') # pdb
	reflines1 = ifile1.readlines()
	reflines2 = ifile2.readlines()

	resid,ffchg = [],[]
	for obj in reflines1:
		line = obj.strip()
		if 'residue' in line or 'Link' in line:
			continue
		else:
			res = line.split()[3]
			atm = line.split()[4]
			rnm = line.split()[2]
			chg = float(line.split()[6])
			resid.append(res+' '+atm+' '+rnm)
			ffchg.append(chg)

	lineN = 0
	for obj in reflines2:
		line  = obj.strip()
		lineN += 1
		if 'POSITION' in line:
			coord1 = lineN
		elif 'BOX' in line:
			coord2 = lineN-1
		else:
			continue

	lineN   = 0
	string  = ''
	for obj in reflines2:
		line  = obj.strip()
		lineN += 1
		if lineN > coord1 and lineN < coord2:
			if 'XXX' in line:
				continue
			elif 'SOL' in line:
				if line.split()[2] == 'OW':
					string  += data_format(line,-0.834)
				else:
					string  += data_format(line,0.417)
			elif ' CL' in line:
				string  += data_format(line,-1.000)
			else:
				res = line.split()[1]+' '+line.split()[2]+' '+line.split()[0]
				chg = ffchg[resid.index(res)]
				string  += data_format(line,chg)
	print(inputfile2,': done')
	ofile = open(inputfile2.split('g96')[0]+'dat','w')
	ofile.write(string)

main(sys.argv[1],sys.argv[2])

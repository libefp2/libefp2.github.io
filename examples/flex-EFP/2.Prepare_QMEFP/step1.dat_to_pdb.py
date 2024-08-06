import sys
'''
python step1.dat_to_pdb.py bchl359-71156.dat
output --> bchl359-71156.pdb
'''
def rewrite(line,lineN):
	x0 = 'ATOM'
	x1 = '{:{align}{width}}'.format('%s'%str(lineN),align='>',width=7)
	x2 = '{:{align}{width}}'.format('%s'%line.split()[1],align='>',width=5)
	if line.split()[2] == 'NVAL':
		resn = 'VAL'
	elif line.split()[2] == 'CGLN':
		resn = 'GLN'
	else:
		resn = line.split()[2]

	x3 = '{:{align}{width}}'.format('%s'%resn,align='>',width=4)
	x4 = '{:{align}{width}}'.format('%s'%line.split()[3],align='>',width=6)
	x5 = '{:{align}{width}}'.format('%.10f'%float(line.split()[4]),align='>',width=18)
	x6 = '{:{align}{width}}'.format('%.10f'%float(line.split()[5]),align='>',width=18)
	x7 = '{:{align}{width}}'.format('%.10f'%float(line.split()[6]),align='>',width=18)
	x8 = '{:{align}{width}}'.format('%s'%'1.00  0.00',align='>',width=12)
	x9 = '{:{align}{width}}'.format('%s'%line.split()[0],align='>',width=12)
	return x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+'\n'

def main(inputfile):
	ifile    = open(inputfile,'r')
	reflines = ifile.readlines()

	lineN  = 0
	string = ''
	for obj in reflines:
		line  = obj.strip()
		lineN += 1
		string += rewrite(line,lineN)

	print(inputfile)
	ofile = open(inputfile.split('dat')[0]+'pdb','w')
	ofile.write(string)

main(sys.argv[1])

import sys,math
'''
python step2.bchl_reassign.py bchl359-71156.dat
output --> new_bchl359-71156.dat
'''
def reassign(line):
	part1,part2,part3,part4 = 19,28,81,139
	repeat                  = 0
	head1,head2,tail1,tail2 = [],[],[],[]

	re_ordered = []
	for i in range(len(line)):
		if i < part1 + repeat:
			head1.append(line[i])
		elif i >= part1 + repeat and i < part2 + repeat:
			tail1.append(line[i])
		elif i >= part2 + repeat and i < part3 + repeat:
			head2.append(line[i])
		else:
			tail2.append(line[i])
		if i == 139 or i == 279 or i == 419 or i == 559 or i == 699 or i == 839 or i == 979 or i == 1119:
			re_ordered              += head1+head2+tail1+tail2
			head1,head2,tail1,tail2  = [],[],[],[]
			repeat                  += 140

	string = ''
	for i in range(len(re_ordered)):
		string += re_ordered[i]+'\n'

	return string


def main(inputfile1):
	ifile           = open(inputfile1,'r')
	reflines        = ifile.readlines()
	bchl,amino,solv = [],[],[]
	temp1,temp2     = [],[]
	for obj in reflines:
		line = obj.strip()
		if line.split()[2] == 'BCL':
			temp1.append(line)
			if 'H203  BCL    366' in line or 'H203  BCL    732' in line or 'H203  BCL   1098' in line:
				bchl.append(temp1)
				temp1 = []
		elif line.split()[2] == 'QSL' or line.split()[2] == 'SOL' or line.split()[2] == 'CL':
			solv.append(line)
		else:
			temp2.append(obj)
			if 'OC2 CGLN    358' in line or 'OC2 CGLN    724' in line or 'OC2 CGLN   1090' in line:
				amino.append(temp2)
				temp2 = []

	bchl_re = []
	for i in range(len(bchl)):
		bchl_re.append(reassign(bchl[i]))


	string = ''
	for i in range(len(amino)):
		for j in range(len(amino[i])):
			string += amino[i][j]
		string += bchl_re[i]

	for i in range(len(solv)-1):
		string += solv[i]+'\n'

	print(inputfile1)
	ofile = open('new_'+inputfile1,'w')
	ofile.write(string + solv[len(solv)-1])

main(sys.argv[1])

import sys,math
'''
python step5.create_mapfile.py new_new_new_bchl359-71156.dat 
output --> bchl359.map
'''
def charge(res):
	negative = ['GLU','ASP']
	positive = ['LYS','ARG','HIP']
	if res in negative:
		return '-1'
	elif res in positive:
		return '1'
	else:
		return '0'

def optimized_aa(pn):
	if pn == 359:
		return [103,142,144]
	elif pn == 360:
		return [65,862,863]
	elif pn == 361:
		return [8,33,290]
	elif pn == 362:
		return [282,286,337]
	elif pn == 363:
		return [47,234,235,236]
	elif pn == 364:
		return [88,138,176]
	elif pn == 365:
		return [190,289,290]
	elif pn == 366:
		return [154,168,848,849,850,851]
	else:
		print('not available site!', str(pn))
		exit()

def write_mapfile(pdb,pn):
	opt    = optimized_aa(pn)
	string = ''
	for i in range(len(pdb)):
		line = pdb[i].strip()
		if line.split()[2] == 'QSL':
			continue
		elif line.split()[8] == 'QM':
			continue
		elif line.split()[8] == 'MM':
			continue
		elif line.split()[8] == 'EFP':
			if line.split()[2] == 'BCL':
				if line.split()[1] == 'MG':
					string += '$residue\n'
					string += 'name = h_'+line.split()[3]+'\n'
					string += 'preatoms = '+str(i+1)+'-'+str(i+1+71)+'\n'
					string += 'ch1 = '+str(i+1+10)+','+str(i+1+72)+'\n'
					string += 'postatoms = '+str(i+1)+'-'+str(i+1+71)+'\n'
					string += 'rescharge = 0\n'
					string += 'usefp = h_'+line.split()[3]+'\n'
					string += '$end\n\n'
					string += '$residue\n'
					string += 'name = t_'+line.split()[3]+'\n'
					string += 'preatoms = '+str(i+1+72)+'-'+str(i+1+139)+'\n'
					string += 'ch1 = '+str(i+1+72)+','+str(i+1+10)+'\n'
					string += 'postatoms = '+str(i+1+72)+'-'+str(i+1+139)+'\n'
					string += 'rescharge = 0\n'
					string += 'usefp = t_'+line.split()[3]+'\n'
					string += '$end\n\n'
				else:
					continue
			elif line.split()[2] == 'NVAL':
				if line.split()[1] == 'N':
					string += '$residue\n'
					string += 'name = v_'+line.split()[3]+'\n'
					string += 'preatoms = '+str(i+1)+'-'+str(i+1+15)+'\n'
					string += 'ch1 = '+str(i+1+4)+','+str(i+1+16)+'\n'
					string += 'postatoms = '+str(i+1)+'-'+str(i+1+15)+'\n'
					string += 'rescharge = 1\n'
					string += 'usefp = v_'+line.split()[3]+'\n'
					string += '$end\n\n'
				else:
					continue
			elif line.split()[2] == 'CGLN':
				if line.split()[1] == 'N':
					string += '$residue\n'
					string += 'name = g_'+line.split()[3]+'\n'
					string += 'preatoms = '+str(i+1-2)+'-'+str(i+1+17)+'\n'
					string += 'ch1 = '+str(i+1-2)+','+str(i+1+16)+'\n'
					string += 'postatoms = '+str(i+1-2)+'-'+str(i+1-8)+'\n'
					string += 'rescharge = -1\n'
					string += 'usefp = g_'+line.split()[3]+'\n'
					string += '$end\n\n'
				else:
					continue
			else:
				if line.split()[1] == 'N':
					for j in range(i,i-50,-1):
						if pdb[j].strip().split()[1] == 'C':
							c1  = j+1
						elif pdb[j].strip().split()[1] == 'CA':
							ca1 = j+1
							break
						else:
							continue
				elif line.split()[1] == 'CA':
					ca2 = i+1
				elif line.split()[1] == 'C':
					c2  = i+1
					if int(line.split()[3]) in opt:
						if int(line.split()[3]) == 234 and pn == 363:
							continue
						elif int(line.split()[3]) == 848 and pn == 366:
							continue
						elif int(line.split()[3]) == 290 and pn == 365:
							string += '$residue\n'
							string += 'name = apt-'+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'
							string += 'preatoms = '+str(c1)+'-'+str(c2-1)+'\n'
							string += 'ch1 = '+str(ca2)+','+str(c2)+'\n'
							string += 'ch2 = '+str(c1)+','+str(ca1)+'\n'
							string += 'postatoms = '+str(c1)+'-'+str(c2-1)+'\n'
							string += 'rescharge = '+charge(line.split()[2])+'\n'
							string += 'usefp = apt-'+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'
							string += '$end\n\n'
						else:
							string += '$residue\n'
							string += 'name = opt-'+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'
							string += 'preatoms = '+str(c1)+'-'+str(c2-1)+'\n'
							string += 'ch1 = '+str(ca2)+','+str(c2)+'\n'
							string += 'ch2 = '+str(c1)+','+str(ca1)+'\n'
							string += 'postatoms = '+str(c1)+'-'+str(c2-1)+'\n'
							string += 'rescharge = '+charge(line.split()[2])+'\n'
							string += 'usefp = opt-'+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'                                   
							string += '$end\n\n'
					else:
						string += '$residue\n'
						string += 'name = '+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'
						string += 'preatoms = '+str(c1)+'-'+str(c2-1)+'\n'
						string += 'ch1 = '+str(ca2)+','+str(c2)+'\n'
						string += 'ch2 = '+str(c1)+','+str(ca1)+'\n'
						string += 'postatoms = '+str(c1)+'-'+str(c2-1)+'\n'
						string += 'rescharge = '+charge(line.split()[2])+'\n'
						string += 'usefp = '+line.split()[2][0:1].lower()+'_'+line.split()[3]+'\n'                                   
						string += '$end\n\n'


				else:
					continue
		else:
			print('no environment?')
			exit()
		
	return string 

def main(inputfile):
	ifile    = open(inputfile,'r')
	reflines = ifile.readlines()
	pn       = int(inputfile.split('-')[0].split('bchl')[1])
	pdb      = []
	for obj in reflines:
		line = obj.strip()
		pdb.append(line)

	mapfile = write_mapfile(pdb,pn)
	ofile = open('bchl'+str(pn)+'.map','w')
	ofile.write(mapfile)

main(sys.argv[1])

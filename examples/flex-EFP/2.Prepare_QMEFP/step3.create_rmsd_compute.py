import sys
import numpy as np

def get_coordinates_inp(filename):
	ifile    = open(filename,'r')
	reflines = ifile.readlines()
	lineN    = 0
	for obj in reflines:
		lineN += 1
		if '$data' in obj.strip():
			cor1 = lineN+2
		elif 'comment Atoms' in obj.strip():
			cor2 = lineN-1
		else:
			continue
	lineN = 0
	atoms,coors = [],[]
	for obj in reflines:
		line  = obj.strip()
		lineN += 1
		temp  = []
		if lineN > cor1 and lineN < cor2:
			atoms.append(line.split()[0][0:1])
			temp.append(float(line.split()[2]))
			temp.append(float(line.split()[3]))
			temp.append(float(line.split()[4]))
			coors.append(temp)
		else:
			continue

	return np.array(atoms),np.array(coors)

def get_coordinates_efp(filename):
	ifile    = open(filename,'r')
	reflines = ifile.readlines()
	conv     = 0.529177249
	lineN    = 0
	for obj in reflines:
		lineN += 1
		if 'COORDINATES (BOHR)' in obj.strip():
			cor1 = lineN
		elif 'BO' in obj.strip().split()[0]:
			cor2 = lineN
			break
		else:
			continue
	lineN = 0
	atoms,coors = [],[]
	for obj in reflines:
		line  = obj.strip()
		lineN += 1
		temp  = []
		if lineN > cor1 and lineN < cor2:
			atoms.append(line.split()[0][3])
			temp.append(float(line.split()[1])*conv)
			temp.append(float(line.split()[2])*conv)
			temp.append(float(line.split()[3])*conv)
			coors.append(temp)
		else:
			continue
	return np.array(atoms),np.array(coors)

def atom_number(A):
	temp = list()
	for i in range(len(A)):
		temp.append(number(A[i]))
	return temp

def number(atom):
	if atom == "H":
		amass = 1.0078250
	elif atom == "C":
		amass = 12.0000000
	elif atom == "N":
		amass = 14.0030700
	elif atom == "O":
		amass = 15.9949100
	elif atom == "S":
		amass = 32.0650000
	else:
		print(atom)
		exit("not standard aminoacids add atomic mass information")
	return amass

def centroid(X):
	"""
	Calculate the centroid from a vectorset X
	"""
	C = sum(X)/len(X)
	return C

def kabsch_rmsd(P,Q,A):
	"""
	Rotate matrix P unto Q and calculate the RMSD
	"""
	P = kabsch_rotate(P,Q)
	return rmsd(P,Q,A)

def kabsch_rotate(P,Q):
	"""
	Rotate matrix P unto matrix Q using Kabsch algorithm
	"""
	U = kabsch(P,Q)
	# Rotate P
	P = np.dot(P,U)
	return P

def kabsch(P,Q):
	"""
	The optimal rotation matrix U is calculated and then used to rotate matrix
	P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
	calculated.
	
	Using the Kabsch algorithm with two sets of paired point P and Q,
	centered around the center-of-mass.
	Each vector set is represented as an NxD matrix, where D is the
	the dimension of the space.
	
	The algorithm works in three steps:
	- a translation of P and Q
	- the computation of a covariance matrix C
	- computation of the optimal rotation matrix U
	
	http://en.wikipedia.org/wiki/Kabsch_algorithm
	
	Parameters:
	P -- (N, number of points)x(D, dimension) matrix
	Q -- (N, number of points)x(D, dimension) matrix
	
	Returns:
	U -- Rotation matrix
	"""
	# Computation of the covariance matrix
	C = np.dot(np.transpose(P), Q)
	
	# Computation of the optimal rotation matrix
	# This can be done using singular value decomposition (SVD)
	# Getting the sign of the det(V)*(W) to decide
	# whether we need to correct our rotation matrix to ensure a
	# right-handed coordinate system.
	# And finally calculating the optimal rotation matrix U
	# see http://en.wikipedia.org/wiki/Kabsch_algorithm
	V, S, W = np.linalg.svd(C)
	d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

	if d:
		S[-1] = -S[-1]
		V[:, -1] = -V[:, -1]
	
	# Create Rotation matrix U
	U = np.dot(V, W)
	
	return U

def rmsd(V,W,A):
	"""
	Calculate Root-mean-square deviation from two sets of vectors V and W.
	"""
	D = len(V[0])
	N = len(V)
	rmsd,mass,j = 0.0,0.0,0
	for v, w in zip(V, W):
		rmsd += A[j]*sum([(v[i]-w[i])**2.0 for i in range(D)])
		mass += A[j]
		j    += 1
	return np.sqrt(rmsd)/np.sqrt(mass)


def main(inputfile1,inputfile2,inputfile3):
	ifile1    = open(inputfile1,'r') # amino list
	ifile2    = open(inputfile2,'r') # input list
	ifile3    = open(inputfile3,'r') # efpdb list
	reflines1 = ifile1.readlines()
	reflines2 = ifile2.readlines()
	reflines3 = ifile3.readlines()
	pig1      = [359,360,361,362,363,364,365,366]
	pig2      = [725,726,727,728,729,730,731,732]
	pig3      = [1091,1092,1093,1094,1095,1096,1097,1098]
	pig       = pig1+pig2+pig3
	resn,resd = [],[]
	for obj in reflines1:
		line = obj.strip()
		resn.append(line.split()[0].lower())
		resd.append(int(line.split()[1]))
	efrg,aass = [],[]
	for obj in reflines2:
		line = obj.strip()
		if int(line.split()[0].split('.')[0].split('_')[-1]) in pig:
			continue
		else:
			efrg.append(line.split()[0])
			aass.append(resn[resd.index(int(line.split()[0].split('.')[0].split('_')[-1]))])
	data = []
	path = '/group/lslipche/data/yb_boss/flexible_efp/efpdb/'
	for obj in reflines3:
		line = obj.strip()
		data.append(line.split()[0])

	if 'brown' in path:
		core = str(24)
	else:
		core = str(20)

	string1,string2 = '',''
	for i in range(len(efrg)):
		aa1           = aass[i]
		p_atoms,p_all = get_coordinates_inp(efrg[i])
		P,A           = p_all,p_atoms
		A             = atom_number(A)
		Pc            = centroid(P)
		P            -= Pc
		temp1,temp2   = [],[]
		for j in range(len(data)):
			if 'nval' in data[j]:
				aa2 = 'nval'
			else:
				aa2 = data[j][0:3]
			if aa1 == aa2:
				q_atoms,q_all = get_coordinates_efp(path+aa1+'/'+data[j])
				if np.count_nonzero(p_atoms != q_atoms):
					exit("Atoms not in the same order")
				Q   = q_all
				Qc  = centroid(Q)
				Q  -= Qc
#				A   = atom_number(A)
				temp1.append(kabsch_rmsd(P,Q,A))
				temp2.append(path+aa1+'/'+data[j])
		string2 += temp2[temp1.index(min(temp1))]+' '+efrg[i]+' '+str(min(temp1))+'\n'
		# THRESHOLD in Angstrom
		if min(temp1) < 0.2:
			string1 += 'python step4.Flexible_V5.py '+temp2[temp1.index(min(temp1))]+' '+efrg[i]+' -d -xr\n'
		else:
			string1 += 'gms_slurm -p '+str(core)+' -q standby -w 00:30:00 -test '+efrg[i]+'\n'
		print(str(i+1)+' / '+ str(len(efrg))+' ....',aa1,temp2[temp1.index(min(temp1))].split('/')[-2])

#	string1 += 'python step4.Flexible_V5-head.py head.efp h_359.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_360.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_361.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_362.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_363.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_364.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_365.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_366.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_731.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_1095.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_1096.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5-head.py head.efp h_1098.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_1095.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_1096.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_1098.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_359.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_360.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_361.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_362.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_363.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_364.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_365.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_366.inp -d -xr\n'
#	string1 += 'python step4.Flexible_V5.py tail.efp t_731.inp -d -xr\n'
	string1 += 'rm -f ~/scr/*.script\n'

	ofile1 = open('step4.Flexible_V5.sh','w')
	ofile2 = open('step3.create_rmsd_compute.dat','w')
	ofile1.write(string1)
	ofile2.write(string2)

main(sys.argv[1],sys.argv[2],sys.argv[3])

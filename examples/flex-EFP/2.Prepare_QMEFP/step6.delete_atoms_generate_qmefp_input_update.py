import sys,math
'''
python step6.gen_qmefp.py .map .dat
'''
def add_hydrogen(xyz1,xyz2):
	x1,y1,z1 = float(xyz1.split()[4]),float(xyz1.split()[5]),float(xyz1.split()[6])
	x2,y2,z2 = float(xyz2.split()[4]),float(xyz2.split()[5]),float(xyz2.split()[6])
	
	if xyz2.split()[0] == 'N':
		desired_length = 1.00339
	else:
		desired_length = 1.07886
	
	actual__length = math.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))
	x3             = ((x1-x2)*desired_length/actual__length)+x2
	y3             = ((y1-y2)*desired_length/actual__length)+y2
	z3             = ((z1-z2)*desired_length/actual__length)+z2
	return pdb_format(x3,y3,z3)

def pdb_format(hx,hy,hz):
	x1 = '{:{align}{width}}'.format('%s'%'H',align='<',width=2)
	x2 = '{:{align}{width}}'.format('%s'%'H000',align='>',width=5)
	x3 = '{:{align}{width}}'.format('%s'%'XXX',align='>',width=5)
	x4 = '{:{align}{width}}'.format('%s'%'XXX',align='>',width=7)
	x5 = '{:{align}{width}}'.format('%.8f'%(float(hx)),align='>',width=15)
	x6 = '{:{align}{width}}'.format('%.8f'%(float(hy)),align='>',width=15)
	x7 = '{:{align}{width}}'.format('%.8f'%(float(hz)),align='>',width=15)
	x8 = '{:{align}{width}}'.format('%.10f'%float(0.0),align='>',width=15)
	x9 = '{:{align}{width}}'.format('%s'%'EFP',align='>',width=5)
	return x1+x2+x3+x4+x5+x6+x7+x8+x9

def qm_format(line):
	string = ''
	for i in range(len(line)):
		x1 = '{:{align}{width}}'.format('%s'%line[i].split()[0],align='<',width=3)
		x2 = '{:{align}{width}}'.format('%s'%figure_atomic_number(line[i].split()[0]),align='>',width=6)
		x3 = '{:{align}{width}}'.format('%.8f'%float(line[i].split()[4]),align='>',width=15)
		x4 = '{:{align}{width}}'.format('%.8f'%float(line[i].split()[5]),align='>',width=15)
		x5 = '{:{align}{width}}'.format('%.8f'%float(line[i].split()[6]),align='>',width=15)
		string += x1+x2+x3+x4+x5+'\n'
	return string+' $end\n'

def compute_distance(x1,y1,z1,x2,y2,z2):
	dx,dy,dz = float(x2)-float(x1),float(y2)-float(y1),float(z2)-float(z1)
	r        = math.sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2))
	return r

def figure_atomic_number(a):
	if a.lower() == 'mg':
		return '12.0'
	elif a.lower() == 'h':
		return '1.0'
	elif a.lower() == 'o':
		return '8.0'
	elif a.lower() == 'n':
		return '7.0'
	elif a.lower() == 'c':
		return '6.0'
	elif a.lower() == 's':
		return '16.0'
	else:
		print('What?')
		exit()

def qm_option(chg):
	string =  ' $contrl scftyp=rhf runtyp=energy ispher=1 dfttyp=pbe0\n'
	string += ' icut=12 maxit=200 tddft=excite exetyp=run icharg=0 $end\n'
	string += ' $system mwords=1000 $end\n'
	string += ' $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n'
	string += ' $guess guess=huckel $end\n'
	string += ' $tddft nstate=4 iroot=1 mult=1 nrad=96 nleb=302 $end\n'
	string += ' $scf diis=.t. fdiff=.f. dirscf=.t. $end\n'
	string += ' $DATA\n'
	string += ' bchl\n'
	string += 'C1\n'
	return string

def assign_qm_atoms_xyz(reflines2):
	qxyz   = []
	for obj in reflines2:
		line = obj.strip()
		if line.split()[8] == 'QM':
			qxyz.append(line)
		else:
			continue

	qm_atoms =  qm_option(0)
	qm_atoms += qm_format(qxyz)

	return qm_atoms

def assign_qm_atoms_add(reflines2,pn):
	hid           = ['CB','HB1','HB2','CG','ND1','HD1','CE1','HE1','NE2','CD2','HD2']
	hie           = ['CB','HB1','HB2','CG','ND1','CE1','HE1','NE2','HE2','CD2','HD2']
	his_h1,his_h2 = ['CA'],['CB']

	pre_pep       = ['C','O']
	pre_h1,pre_h2 = ['CA'],['C']
	leu           = ['N','H','CA','HA','CB','HB1','HB2','CG','HG','CD1','HD11','HD12','HD13','CD2','HD21','HD22','HD23','C','O']
	tyr           = ['N','H','CA','HA','CB','HB1','HB2','CG','CD1','HD1','CE1','HE1','CZ','OH','HH','CE2','HE2','CD2','HD2','C','O']
	nex_pep       = ['N','H','CA','HA']
	nex_h1,nex_h2 = ['CB','C'],['CA','CA']

	if pn == 359:
		add_to_aa  = {'103':hid}
		add_h000_1 = {'103':his_h1}
		add_h000_2 = {'103':his_h2}
	elif pn == 361:
		add_to_aa  = {'290':hid}
		add_h000_1 = {'290':his_h1}
		add_h000_2 = {'290':his_h2}
	elif pn == 362:
		add_to_aa  = {'282':hid}
		add_h000_1 = {'282':his_h1}
		add_h000_2 = {'282':his_h2}
	elif pn == 363:
		add_to_aa  = {'233':pre_pep,'234':leu,'235':nex_pep}
		add_h000_1 = {'233':pre_h1,'234':[],'235':nex_h1}
		add_h000_2 = {'233':pre_h2,'234':[],'235':nex_h2}
	elif pn == 364:
		add_to_aa  = {'138':hid}
		add_h000_1 = {'138':his_h1}
		add_h000_2 = {'138':his_h2}
	elif pn == 365:
		add_to_aa  = {'289':hie}
		add_h000_1 = {'289':his_h1}
		add_h000_2 = {'289':his_h2}
	elif pn == 366:
		add_to_aa  = {'115':pre_pep,'116':tyr,'117':nex_pep}
		add_h000_1 = {'115':pre_h1,'116':[],'117':nex_h1}
		add_h000_2 = {'115':pre_h2,'116':[],'117':nex_h2}
	elif pn == 1098:
		add_to_aa  = {'847':pre_pep,'848':tyr,'849':nex_pep}
		add_h000_1 = {'847':pre_h1,'848':[],'849':nex_h1}
		add_h000_2 = {'847':pre_h2,'848':[],'849':nex_h2}
	else:
		print('not possible pigment')
		exit()

	hxyz1,hxyz2,axyz,bxyz = [],[],[],[]
	for obj in reflines2:
		line = obj.strip()
		if 'SOL' in line or 'OUT' in line or 'INS' in line or 'QSL' in line:
			continue
		elif line.split()[2] == 'BCL':
			if int(line.split()[3]) == pn:
				bxyz.append(line)
			else:
				continue
		else:
			if line.split()[3] in add_to_aa.keys():
				if line.split()[1] in add_h000_1[line.split()[3]]:
					hxyz1.append(line)
				elif line.split()[1] in add_h000_2[line.split()[3]]:
					hxyz2.append(line)
					axyz.append(line)
				elif line.split()[1] in add_to_aa[line.split()[3]]:
					axyz.append(line)
				else:
					continue
			else:
				continue

	if len(hxyz1) == len(hxyz2):
		for i in range(len(hxyz1)):
			axyz.append(add_hydrogen(hxyz1[i],hxyz2[i]))
	elif len(hxyz1) > len(hxyz2):
		for i in range(len(hxyz2)):
			axyz.append(add_hydrogen(hxyz1[i],hxyz2[i]))
		axyz.append(add_hydrogen(hxyz1[-1],hxyz2[-1]))
	else:
		for i in range(len(hxyz1)):
			axyz.append(add_hydrogen(hxyz1[i],hxyz2[i]))
		axyz.append(add_hydrogen(hxyz1[-1],hxyz2[-1]))

	qm_atoms  = qm_option(0)
	qm_atoms += qm_format(axyz+bxyz)

	return qm_atoms

def rm_to_atoms(fragname,efparm,coord1,coord2,get_rid_of_aa,buff,rm_to_line):
	h000,hx,hy,hz = [],[],[],[]
	temp1,temp2   = [],[]
	temp3,temp4   = [],[]
	junk1,buf0    = [],[]
	if fragname == 'opt-h_138':
	#	print(fragname)
		for i in range(coord1,coord2-1,1):
			line = efparm[i].strip()
			if 'A05C' in line or 'A07C' in line:
				junk1.append(line)
			else:
				continue
		dx = (float(junk1[0].split()[1])+float(junk1[1].split()[1]))/2.0
		dy = (float(junk1[0].split()[2])+float(junk1[1].split()[2]))/2.0
		dz = (float(junk1[0].split()[3])+float(junk1[1].split()[3]))/2.0
		hx.append(dx);hy.append(dy);hz.append(dz)

	for i in range(coord1,coord2-1,1):
		line = efparm[i].strip()
		if line.split()[0][0:4] in get_rid_of_aa:
			h000.append(line.split()[0])
			temp1.append(int(line.split()[0][1:3]))
			rm_to_line.append(i)
			if 'N' in line.split()[0] or 'O' in line.split()[0] or 'S' in line.split()[0]:
				hx.append(float(line.split()[1]))
				hy.append(float(line.split()[2]))
				hz.append(float(line.split()[3]))
			else:
				continue
		elif line.split()[0][0:4] in buff:
			buf0.append(line.split()[0])
			hx.append(float(line.split()[1]))
			hy.append(float(line.split()[2]))
			hz.append(float(line.split()[3]))
			temp2.append(int(line.split()[0][1:3]))
		elif 'H000' in line.split()[0]:
			h000.append(line.split()[0])
			temp1.append(int(line.split()[0][1:3]))
			rm_to_line.append(i)
			temp3.append('BO'+str(int(efparm[i].strip().split()[0].split('H000')[0].split('A')[1])))
		elif 'BO' in line.split()[0]:
			bmid = line.split()[0].split('BO')[1]
			if len(bmid) == 2:
				atom1,atom2 = int(bmid[0:1]),int(bmid[1:2])
			elif len(bmid) == 3:
				atom1,atom2 = int(bmid[0:2]),int(bmid[2:3])
			else:
				if int(bmid[0:2]) > int(bmid[2:4]):
					atom1,atom2 = int(bmid[0:2]),int(bmid[2:4])
				else:
					atom1,atom2 = int(bmid[0:3]),int(bmid[3:4])
			if atom1 in temp1 or atom2 in temp1:
				h000.append(line.split()[0])
				temp4.append(line.split()[0])
				hx.append(float(line.split()[1]))
				hy.append(float(line.split()[2]))
				hz.append(float(line.split()[3]))
				rm_to_line.append(i)
			elif atom1 in temp2 or atom2 in temp2:
				buf0.append(line.split()[0])
				hx.append(float(line.split()[1]))
				hy.append(float(line.split()[2]))
				hz.append(float(line.split()[3]))
			else:
				continue
	return h000,buf0,hx,hy,hz,rm_to_line,temp3,temp4

def charge_saturation(line,charge):
	neut = float(line.split()[1])+charge
	x0   = '{:{align}{width}}'.format('%s'%line.strip().split()[0],align='<',width=len(line.strip().split()[0]))
	x1   = '{:{align}{width}}'.format('%.10f'%float(neut),align='>',width=23-len(line.strip().split()[0]))
	x2   = '{:{align}{width}}'.format('%.5f'%float(line.strip().split()[2]),align='>',width=10)
	return x0+x1+x2+'\n'

def rm_to_mopole(efparm,h000,saturate,monop1,monop2,rm_to_line):
	temp1,temp2  = [],[]
	for i in range(monop1,monop2-1,1):
		line = efparm[i].strip()
		if line.split()[0] in h000:
			if 'BO' in line.split()[0]:
				rm_to_line.append(i)
				temp2.append(float(line.split()[1])+float(line.split()[2]))
			else:
				rm_to_line.append(i)
				temp1.append(float(line.split()[1])+float(line.split()[2]))
		else:
			continue

	c_chg,n_chg = 0.0,0.0

	if 'BO' not in saturate[0]:
		if len(saturate) == 1:
			for i in range(len(temp1)):
				c_chg += temp1[i]
			for i in range(len(temp2)):
				c_chg += temp2[i]
			chg = [c_chg]
			for i in range(len(saturate)):
				efparm[monop1+int(saturate[i][1:3])-1] = charge_saturation(efparm[monop1+int(saturate[i][1:3])-1],chg[i])
		elif len(saturate) == 2:
			c = [-1]
			for i in range(len(temp1)):
				ind = (len(temp1)-i-1)-len(temp1)
				if ind in c:
					c_chg += temp1[(len(temp1)-i-1)]
				else:
					n_chg += temp1[(len(temp1)-i-1)]
			for i in range(len(temp2)):
				ind = (len(temp2)-i-1)-len(temp2)
				if ind in c:
					c_chg += temp2[(len(temp2)-i-1)]
				else:
					n_chg += temp2[(len(temp2)-i-1)]
			chg = [c_chg,n_chg]
			for i in range(len(saturate)):
				efparm[monop1+int(saturate[i][1:3])-1] = charge_saturation(efparm[monop1+int(saturate[i][1:3])-1],chg[i])
		else:
			c = [-1]
			for i in range(len(temp1)):
				ind = (len(temp1)-i-1)-len(temp1)
				if ind in c:
					c_chg += temp1[(len(temp1)-i-1)]
				else:
					n_chg += temp1[(len(temp1)-i-1)]
			for i in range(len(temp2)):
				ind = (len(temp2)-i-1)-len(temp2)
				if ind in c:
					c_chg += temp2[(len(temp2)-i-1)]
				else:
					n_chg += temp2[(len(temp2)-i-1)]
			chg = [c_chg,n_chg/2.0,n_chg/2.0]
			for i in range(len(saturate)):
				efparm[monop1+int(saturate[i][1:3])-1] = charge_saturation(efparm[monop1+int(saturate[i][1:3])-1],chg[i])
	else:
		if len(h000) == 2:
			efparm[monop1+(int(h000[1].split(saturate[0])[1])-1)] = charge_saturation(efparm[monop1+int(h000[1].split(saturate[0])[1])-1],temp1[0]+temp2[0])
		elif len(h000) == 4:
			efparm[monop1+(int(h000[2].split(saturate[0])[1])-1)] = charge_saturation(efparm[monop1+int(h000[2].split(saturate[0])[1])-1],temp1[0]+temp2[0])
			efparm[monop1+(int(h000[3].split(saturate[1])[1])-1)] = charge_saturation(efparm[monop1+int(h000[3].split(saturate[1])[1])-1],temp1[1]+temp2[1])
		elif len(h000) == 6:
			efparm[monop1+(int(h000[3].split(saturate[0])[1])-1)] = charge_saturation(efparm[monop1+int(h000[3].split(saturate[0])[1])-1],temp1[0]+temp2[0])
			efparm[monop1+(int(h000[4].split(saturate[1])[1])-1)] = charge_saturation(efparm[monop1+int(h000[4].split(saturate[1])[1])-1],temp1[1]+temp2[1])
			efparm[monop1+(int(h000[5].split(saturate[2])[1])-1)] = charge_saturation(efparm[monop1+int(h000[5].split(saturate[2])[1])-1],temp1[2]+temp2[2])
		else:
			print('more than three additional hydrogen??')
			exit()
	
	tot = 0.0
	for i in range(monop1,monop2-1,1):
		if efparm[i].strip().split()[0] in h000:
			continue
		else:
			tot += (float(efparm[i].strip().split()[1])+float(efparm[i].strip().split()[2]))

	return rm_to_line,efparm,tot

def rm_to_dipole(efparm,h000,dipol1,dipol2,rm_to_line):
	for i in range(dipol1,dipol2-1,1):
		if efparm[i].strip().split()[0] in h000:
			rm_to_line.append(i)
		else:
			continue
	return rm_to_line

def rm_to_qupole(efparm,h000,qupol1,qupol2,rm_to_line):
	for i in range(qupol1,qupol2-1,1):
		if efparm[i].strip().split()[0] in h000:
			rm_to_line.append(i)
			rm_to_line.append(i+1)
		else:
			continue
	return rm_to_line

def rm_to_ocpole(efparm,h000,ocpol1,ocpol2,rm_to_line):
	for i in range(ocpol1,ocpol2-1,1):
		if efparm[i].strip().split()[0] in h000:
			rm_to_line.append(i)
			rm_to_line.append(i+1)
			rm_to_line.append(i+2)
		else:
			continue
	return rm_to_line

def rm_to_static(efparm,hx,hy,hz,statc1,statc2,rm_to_line):
	temp1,temp2 = [],[]
	for i in range(statc1,statc2-1,1):
		if 'CT' in efparm[i].strip().split()[0]:
			temp1.append(efparm[i].strip())
		else:
			continue
	
	for i in range(len(hx)):
		junk = []
		for j in range(len(temp1)):
			dx = hx[i] - float(temp1[j].split()[1])
			dy = hy[i] - float(temp1[j].split()[2])
			dz = hz[i] - float(temp1[j].split()[3])
			r2 = dx**2+dy**2+dz**2
			junk.append(r2)
		temp2.append(temp1[junk.index(min(junk))].split()[0])
		min_r1 = min(junk)
		junk[junk.index(min(junk))] = 10000
		min_r2 = min(junk)
		if min_r2- min_r1 <= 0.5:
			temp2.append(temp1[junk.index(min(junk))].split()[0])
		else:
			continue

	for i in range(statc1,statc2-1,1):
		if efparm[i].strip().split()[0] in temp2:
			rm_to_line.append(i)
			rm_to_line.append(i+1)
			rm_to_line.append(i+2)
			rm_to_line.append(i+3)
		else:
			continue
	return rm_to_line,temp2

def rm_to_screen(efparm,h000,screen,end,rm_to_line):
	for i in range(screen,end-1,1):
		if efparm[i].strip().split()[0] in h000:
			rm_to_line.append(i)
		else:
			continue
	return rm_to_line

def h_coord(reflines,fragname,get_rid_of_aa,buff,saturate):
	efparm = []
	lineN  = 0
	for obj in reflines:
		line  = obj.strip()
		lineN += 1
		efparm.append(obj)
		if 'COORDINATES (BOHR)' in line:
			coord1 = lineN
		elif 'MONOPOLES' in line:
			coord2 = lineN-1
			monop1 = lineN
		elif 'DIPOLES' in line:
			monop2 = lineN-1
			dipol1 = lineN
		elif 'QUADRUPOLES' in line:
			dipol2 = lineN-1
			qupol1 = lineN
		elif 'OCTUPOLES' in line:
			qupol2 = lineN-1
			ocpol1 = lineN
		elif 'POLARIZABLE POINTS' in line:
			ocpol2 = lineN-1
			statc1 = lineN
		elif 'SCREEN2' in line:
			statc2 = lineN-1
			screen = lineN
		elif line.split()[0] == 'SCREEN':
			screen2 = lineN-1
			scr     = lineN
		else:
			continue

	rm_to_line = []
	h000,buf0,hx,hy,hz,rm_to_line,connect,junk = rm_to_atoms(fragname,efparm,coord1,coord2,get_rid_of_aa,buff,rm_to_line)
	if len(saturate) > 0:
		rm_to_line,efparm,monop_chg   = rm_to_mopole(efparm,h000,saturate,monop1,monop2,rm_to_line)
	else:
		rm_to_line,efparm,monop_chg   = rm_to_mopole(efparm,h000,connect,monop1,monop2,rm_to_line)
	rm_to_line                            = rm_to_dipole(efparm,h000+buf0,dipol1,dipol2,rm_to_line)
	rm_to_line                            = rm_to_qupole(efparm,h000+buf0,qupol1,qupol2,rm_to_line)
	rm_to_line                            = rm_to_ocpole(efparm,h000+buf0,ocpol1,ocpol2,rm_to_line)
	rm_to_line,rm_to_lmo                  = rm_to_static(efparm,hx,hy,hz,statc1,statc2,rm_to_line)
	rm_to_line                            = rm_to_screen(efparm,h000,screen,lineN,rm_to_line)

	parm_removed = []

	parm_removed = []
#	if len(saturate) > 0:
#		wrt = scr-1
#	else:
#		wrt = len(efparm)-1
	wrt = scr-1

	for i in range(wrt):
		line = efparm[i].strip()
		if i in rm_to_line:
			continue
		else:
			if '$FRAGNAME' in line:
				parm_removed.append(' $'+fragname+'\n')
			else:
				parm_removed.append(efparm[i])
	parm_removed.append(' $END\n')
	return parm_removed,monop_chg

def assign_ef_atoms(reflines2,frag,parm,pn):
	temp1,temp2 = list(),list()
	water       = ['A01O1','A02H2','A03H3']
	total       = 0.0

	hid         = ['A05C','A07C','A08H','A09H','A10C','A11N','A12H','A13C','A14H','A15N','A16C','A17H','A18H','A19H']
	hie         = ['A05C','A07C','A08H','A09H','A10C','A11N','A12C','A13H','A14N','A15H','A16C','A17H','A18H','A19H']
	his_buf     = ['A03N','A06H'] # just for lmo deletion
	his_sat     = ['A01C']
	his_nex     = []
	his_nex_buf = ['A01C'] # just for lmo deletion
	his_nex_sat = ['A01C','A05C']
	
	glu         = ['A05C','A16H','A17H']
	glu_buf     = ['A03N','A06H','A07C']
	glu_sat     = ['A01C','A07C']
	phe         = ['A01C','A02O','A03N','A04H','A05C','A06H','A07C','A21H','A22H']
	phe_buf     = ['A08H','A09H','A10C']
	phe_sat     = ['A10C']
	pro         = ['A01C','A15H','A16H']
	pro_buf     = ['A02O','A03N']
	pro_sat     = ['A03N','A13C']
	
	phe2        = ['A05C','A21H','A22H']
	phe2_buf    = ['A03N','A06H','A07C']
	phe2_sat    = ['A01C','A07C']
	tyr2        = ['A01C','A02O','A03N','A04H','A05C','A06H','A07C','A22H','A23H']
	tyr2_buf    = ['A08H','A09H','A10C']
	tyr2_sat    = ['A10C']
	arg         = ['A01C','A25H','A26H']
	arg_buf     = ['A02O','A03N']
	arg_sat     = ['A03N','A05C']

	if pn == 359:
		rm_to_remove = {'opt-h_103':hid,'t_104':his_nex}
		buff_atoms   = {'opt-h_103':his_buf,'t_104':his_nex_buf}
		charge_sat   = {'opt-h_103':his_sat,'t_104':his_nex_sat}
	elif pn == 360:
		rm_to_remove,buff_atoms,charge_sat = {},{},{}
	elif pn == 361:
		rm_to_remove = {'opt-h_290':hid,'g_291':his_nex}
		buff_atoms   = {'opt-h_290':his_buf,'g_291':his_nex_buf}
		charge_sat   = {'opt-h_290':his_sat,'g_291':his_nex_sat}
	elif pn == 362:
		rm_to_remove = {'opt-h_282':hid,'p_283':his_nex}
		buff_atoms   = {'opt-h_282':his_buf,'p_283':his_nex_buf}
		charge_sat   = {'opt-h_282':his_sat,'p_283':['A01C','A13C']}
	elif pn == 363:
		rm_to_remove = {'g_233':glu,'opt-p_235':phe,'opt-p_236':pro}
		buff_atoms   = {'g_233':glu_buf,'opt-p_235':phe_buf,'opt-p_236':pro_buf}
		charge_sat   = {'g_233':glu_sat,'opt-p_235':phe_sat,'opt-p_236':pro_sat}
	elif pn == 364:
		rm_to_remove = {'opt-h_138':hid,'a_139':his_nex}
		buff_atoms   = {'opt-h_138':his_buf,'a_139':his_nex_buf}
		charge_sat   = {'opt-h_138':his_sat,'a_139':his_nex_sat}
	elif pn == 365:
		rm_to_remove = {'opt-h_289':hie,'apt-h_290':his_nex}
		buff_atoms   = {'opt-h_289':his_buf,'apt-h_290':his_nex_buf}
		charge_sat   = {'opt-h_289':his_sat,'apt-h_290':his_nex_sat}
	elif pn == 366:
		rm_to_remove = {'p_115':phe2,'opt-t_117':tyr2,'opt-a_118':arg}
		buff_atoms   = {'p_115':phe2_buf,'opt-t_117':tyr2_buf,'opt-a_118':arg_buf}
		charge_sat   = {'p_115':phe2_sat,'opt-t_117':tyr2_sat,'opt-a_118':arg_sat}
	elif pn == 1098:
		rm_to_remove = {'p_847':phe2,'opt-t_849':tyr2,'opt-a_850':arg}
		buff_atoms   = {'p_847':phe2_buf,'opt-t_849':tyr2_buf,'opt-a_850':arg_buf}
		charge_sat   = {'p_847':phe2_sat,'opt-t_849':tyr2_sat,'opt-a_850':arg_sat}
	else:
		print('not possible pigment')
		exit()

	temp1.append(' $efrag\n')
	temp1.append(' iscrelec=0 iscrpol=1 nodisp nochtr noexrep\n')
	for i in range(len(frag)):
		if int(frag[i].split('_')[1]) == pn:
			continue
		else:
			ifile    = open(frag[i]+'.efp','r')
			reflines = ifile.readlines()
			if frag[i] in rm_to_remove.keys():
				parm_removed,tot_chg = h_coord(reflines,frag[i],rm_to_remove[frag[i]],buff_atoms[frag[i]],charge_sat[frag[i]])
			else:
				parm_removed,tot_chg = h_coord(reflines,frag[i],[],[],[])
			total += tot_chg
			temp1.append('FRAGNAME='+frag[i]+'\n')

			lineN,name,numbering = 0,list(),list()
			for j in range(len(parm_removed)):
				lineN += 1
				if lineN >=6 and lineN <= 8:
					name.append(parm_removed[j].split()[0])
					numbering.append(int(parm_removed[j].split()[0][1:3])-1)
				else:
					continue

			lineN,nummy,x,y,z = 1,0,list(),list(),list()
			for obj in reflines2:
				line  = obj.strip().split()
				if lineN == int(parm[i])+numbering[nummy]:
					x.append(float(line[4]))
					y.append(float(line[5]))
					z.append(float(line[6]))
					nummy += 1
					if nummy == 3:
						nummy = 0
				lineN += 1

			for j in range(len(name)):
				x1 = '{:{align}{width}}'.format('%s'%name[j],align='<',width=8)
				x2 = '{:{align}{width}}'.format('%.6f'%x[j],align='>',width=13)
				x3 = '{:{align}{width}}'.format('%.6f'%y[j],align='>',width=13)
				x4 = '{:{align}{width}}'.format('%.6f'%z[j],align='>',width=13)
				temp1.append(x1+x2+x3+x4+'\n')

			temp2.append(parm_removed)
	
	wx,wy,wz = list(),list(),list()
	for obj in reflines2:
		line = obj.strip().split()
		if line[2] == 'QSL' or line[2] == 'INS':
			if line[8] == 'EFP':
				wx.append(float(line[4]))
				wy.append(float(line[5]))
				wz.append(float(line[6]))
			else:
				continue
		else:
			continue

	count = 0
	for j in range(len(wx)):
		if count == 0:
			temp1.append('FRAGNAME=water\n')
		x1 = '{:{align}{width}}'.format('%s'%water[count],align='<',width=8)
		x2 = '{:{align}{width}}'.format('%.6f'%wx[j],align='>',width=13)
		x3 = '{:{align}{width}}'.format('%.6f'%wy[j],align='>',width=13)
		x4 = '{:{align}{width}}'.format('%.6f'%wz[j],align='>',width=13)
		temp1.append(x1+x2+x3+x4+'\n')
		count += 1
		if count == 3:
			count = 0
	
	water_parm = list()
	ifile1 = open('water.efp','r')
	reflines1 = ifile1.readlines()
	for obj in reflines1:
		line = obj.rstrip()
		water_parm.append(obj)
	
	temp2.append(water_parm)
	return temp1,temp2,total

def assign_mm_atoms(pdb,resn,pn):
	temp1,temp2  = list(),list()
	conv         = 1/0.52917721092
	lineN,neut   = 0,0.0
	name,x,y,z,c = list(),list(),list(),list(),list()
	for i in range(len(pdb)):
		line = pdb[i]
		if line.split()[8] == 'EFP':
			continue
		elif line.split()[8] == 'QM':
			continue
		else:
			if line.split()[2] == 'BCL':
				name.append('O'+str(lineN))
				x.append(float(line.split()[4]))
				y.append(float(line.split()[5]))
				z.append(float(line.split()[6]))
				c.append(float(line.split()[7]))
				lineN += 1
			elif line.split()[2] == 'CL':
				name.append('O'+str(lineN))
				x.append(float(line.split()[4]))
				y.append(float(line.split()[5]))
				z.append(float(line.split()[6]))
				c.append(float(line.split()[7]))
				lineN += 1
			elif line.split()[3] == 'SOL':
				name.append('O'+str(lineN))
				x.append(float(line.split()[4]))
				y.append(float(line.split()[5]))
				z.append(float(line.split()[6]))
				c.append(float(line.split()[7]))
				lineN += 1
			elif line.split()[3] == 'QSL':
				name.append('O'+str(lineN))
				x.append(float(line.split()[4]))
				y.append(float(line.split()[5]))
				z.append(float(line.split()[6]))
				c.append(float(line.split()[7]))
				lineN += 1
			else:
				if line.split()[2] == 'CGLN':
					if int(line.split()[3])-1 in resn:
						if line.split()[1] == 'N':
							name.append('O'+str(lineN))
							x.append(float(pdb[i-2].split()[4]))
							y.append(float(pdb[i-2].split()[5]))
							z.append(float(pdb[i-2].split()[6]))
							c.append(float(pdb[i-2].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(pdb[i-1].split()[4]))
							y.append(float(pdb[i-1].split()[5]))
							z.append(float(pdb[i-1].split()[6]))
							c.append(float(pdb[i-1].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							neut = float(pdb[i-2].split()[7])+float(pdb[i-1].split()[7])
							c.append(float(line.split()[7])-neut)
							lineN += 1
						else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
					else:
						name.append('O'+str(lineN))
						x.append(float(line.split()[4]))
						y.append(float(line.split()[5]))
						z.append(float(line.split()[6]))
						c.append(float(line.split()[7]))
						lineN += 1
				elif line.split()[2] == 'NVAL':
					if int(line.split()[3])+1 in resn:
						if line.split()[1] == 'CA':
							cax = float(line.split()[4])
							cay = float(line.split()[5])
							caz = float(line.split()[6])
							cac = float(line.split()[7])
						elif line.split()[1] == 'C':
							cac += float(line.split()[7])
						elif line.split()[1] == 'O':
							cac += float(line.split()[7])
							name.append('O'+str(lineN))
							x.append(float(cax))
							y.append(float(cay))
							z.append(float(caz))
							c.append(float(cac))
							lineN += 1
						else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
					else:
						name.append('O'+str(lineN))
						x.append(float(line.split()[4]))
						y.append(float(line.split()[5]))
						z.append(float(line.split()[6]))
						c.append(float(line.split()[7]))
						lineN += 1
				else:
					if int(line.split()[3])-1 in resn and int(line.split()[3])+1 not in resn:
						if line.split()[1] == 'N':
							name.append('O'+str(lineN))
							x.append(float(pdb[i-2].split()[4]))
							y.append(float(pdb[i-2].split()[5]))
							z.append(float(pdb[i-2].split()[6]))
							c.append(float(pdb[i-2].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(pdb[i-1].split()[4]))
							y.append(float(pdb[i-1].split()[5]))
							z.append(float(pdb[i-1].split()[6]))
							c.append(float(pdb[i-1].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							neut = float(pdb[i-2].split()[7])+float(pdb[i-1].split()[7])
							c.append(float(line.split()[7])-neut)
							lineN += 1
						else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
					elif int(line.split()[3])-1 not in resn and int(line.split()[3])+1 in resn:
						if line.split()[1] == 'CA':
							cax = float(line.split()[4])
							cay = float(line.split()[5])
							caz = float(line.split()[6])
							cac = float(line.split()[7])
						elif line.split()[1] == 'C':
							cac += float(line.split()[7])
						elif line.split()[1] == 'O':
							cac += float(line.split()[7])
							name.append('O'+str(lineN))
							x.append(float(cax))
							y.append(float(cay))
							z.append(float(caz))
							c.append(float(cac))
							lineN += 1
						else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
					elif int(line.split()[3])-1 in resn and int(line.split()[3])+1 in resn:
						if line.split()[1] == 'N':
							name.append('O'+str(lineN))
							x.append(float(pdb[i-2].split()[4]))
							y.append(float(pdb[i-2].split()[5]))
							z.append(float(pdb[i-2].split()[6]))
							c.append(float(pdb[i-2].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(pdb[i-1].split()[4]))
							y.append(float(pdb[i-1].split()[5]))
							z.append(float(pdb[i-1].split()[6]))
							c.append(float(pdb[i-1].split()[7]))
							lineN += 1
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							neut = float(pdb[i-2].split()[7])+float(pdb[i-1].split()[7])
							c.append(float(line.split()[7])-neut)
							lineN += 1
						elif line.split()[1] == 'CA':
							cax = float(line.split()[4])
							cay = float(line.split()[5])
							caz = float(line.split()[6])
							cac = float(line.split()[7])
						elif line.split()[1] == 'C':
							cac += float(line.split()[7])
						elif line.split()[1] == 'O':
							cac += float(line.split()[7])
							name.append('O'+str(lineN))
							x.append(float(cax))
							y.append(float(cay))
							z.append(float(caz))
							c.append(float(cac))
							lineN += 1
						else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
					else:
							name.append('O'+str(lineN))
							x.append(float(line.split()[4]))
							y.append(float(line.split()[5]))
							z.append(float(line.split()[6]))
							c.append(float(line.split()[7]))
							lineN += 1
	temp1.append('FRAGNAME=PROT')
	for i in range(3):
		x1 = '{:{align}{width}}'.format('%s'%name[i],align='<',width=8)
		x2 = '{:{align}{width}}'.format('%.8f'%x[i],align='>',width=15)
		x3 = '{:{align}{width}}'.format('%.8f'%y[i],align='>',width=15)
		x4 = '{:{align}{width}}'.format('%.8f'%z[i],align='>',width=15)
		temp1.append(x1+x2+x3+x4)
	temp1.append(' $END')
	temp2.append(' $PROT')
	temp2.append('TITLE')
	temp2.append('  COORDINATES (BOHR)')
	for i in range(len(name)):
		x1 = '{:{align}{width}}'.format('%s'%name[i],align='<',width=7)
		x2 = '{:{align}{width}}'.format('%.12f'%float(x[i]*conv),align='>',width=20)
		x3 = '{:{align}{width}}'.format('%.12f'%float(y[i]*conv),align='>',width=20)
		x4 = '{:{align}{width}}'.format('%.12f'%float(z[i]*conv),align='>',width=20)
		temp2.append(x1+x2+x3+x4)
	temp2.append('STOP')
	temp2.append('MONOPOLES')
	for i in range(len(name)):
		x1 = '{:{align}{width}}'.format('%s'%name[i],align='<',width=7)
		x2 = '{:{align}{width}}'.format('%.10f'%float(c[i]),align='>',width=20)
		temp2.append(x1+x2)
	temp2.append('STOP')
	temp2.append('SCREEN2      (FROM VDWSCL=   0.700)')
	for i in range(len(name)):
		x1 = '{:{align}{width}}'.format('%s'%name[i],align='<',width=7)
		x2 = '{:{align}{width}}'.format('%.10f'%float(1.0),align='>',width=15)
		x3 = '{:{align}{width}}'.format('%.10f'%float(10.0),align='>',width=15)
		temp2.append(x1+x2+x3)
#	temp2.append('POLAB 0.3')
	temp2.append('STOP')
	temp2.append(' $END')
	return temp1,temp2,sum(c)
'''
python step6.gen_qmefp.py combined.map _.dat
'''
def main(inputfile1,inputfile2):
	ifile1    = open(inputfile1,'r') # mapfile
	ifile2    = open(inputfile2,'r') # xyz information
	reflines1 = ifile1.readlines()
	reflines2 = ifile2.readlines()
	pn        = int(inputfile2.split('-')[0].split('bchl')[1])
	if pn == 366:
		pn = 1098
	else:
		pn = pn

	frag,resn,parm = list(),list(),list()
	for obj in reflines1:
		line = obj.strip()
		if 'name' in line:
			frag.append(line.split()[2])
			resn.append(int(line.split()[2].split('_')[1]))
		elif 'preatoms' in line:
			parm.append(line.split()[2].split('-')[0])
	pdb = []
	for obj in reflines2:
		line = obj.strip()
		pdb.append(line)

	if pn == 360:
		qm_atoms         = assign_qm_atoms_xyz(reflines2)
	else:
		qm_atoms         = assign_qm_atoms_add(reflines2,pn)

	ef_atoms,ef_parms,ef_chg = assign_ef_atoms(reflines2,frag,parm,pn)
	mm_atoms,mm_parms,mm_chg = assign_mm_atoms(pdb,resn,pn)

	x0 = '{:{align}{width}}'.format('%.5f'%float(ef_chg+mm_chg),align='>',width=8)
	print(inputfile2+' Total FMO charge: '+x0)
	print('-------------')

	ofile = open(inputfile2.split('dat')[0]+'inp','w')
	ofile.write(qm_atoms)
	for i in range(len(ef_atoms)):
		ofile.write(ef_atoms[i])
	for i in range(len(mm_atoms)):
		ofile.write(mm_atoms[i]+'\n')
	for i in range(len(ef_parms)):
		for j in range(len(ef_parms[i])):
			ofile.write(ef_parms[i][j])
	for i in range(len(mm_parms)):
		ofile.write(mm_parms[i]+'\n')
	exit()
#	mm_atoms,mm_param,mm_chg    = assign_mm_atoms(pdb,resn,pn)

main(sys.argv[1],sys.argv[2])

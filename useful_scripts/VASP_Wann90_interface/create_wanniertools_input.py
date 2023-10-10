# Author: Ilias Samathrakis

# The following script creates the input file of wanniertools software based on VASP and Wannier90 output

# Wannier tools: available at 'https://github.com/quanshengwu/wannier_tools'
# VASP: available at 'https://www.vasp.at/'
# Wannier90: available at 'https://wannier.org/'

# The Process requires 
# 1. A self consistent calculation with VASP
# 2. Wannier functions construction with Wannier90
# 3. Band structure calculation with VASP

# DIRECTORIES
# -------------------------------------------------------------------------
wannier90win_dir = ''    # Path of wannier90.win file
wannier90wout_dir = ''   # Path of wannier90.wout file
wanniertools_dir = ''    # Path of wanniertools wt.x file
OUTCAR_dir = ''          # Path of OUTCAR from self consistent calculation
kpts_band_dir = ''       # Path of KPOINTS from band structure calculation
# -------------------------------------------------------------------------

# SUBMISSION FILE
# -------------------------------------------------------------------------
job_name = ''            # Name of job
cores =                  # Number of cores (always integer)
INPUT_dir = ''           # Path of POSCAR
project = ''             # Name of project to submit job
memory =                 # Memory per CPU (2048, 4096, etc)
architecture = ''        # Architecture of HPC system (avx, avx2, avx512, etc)
time = '24:00:00'        # Time (follow this format)
# -------------------------------------------------------------------------

# KPOITS density in wanniertools input file (in respect to KPOINTS of VASP)
# -------------------------------------------------------------------------
density = 1500
# -------------------------------------------------------------------------

# Modify lines 258-269 according to your HPC needs
# Modify lines 271-371 based on the desired wanniertools calculation

import os
import math
import numpy as np
from numpy.linalg import inv

## input file of wanniertools. DO NOT MODIFY
OUTPUT_file = 'wt.in'

lat1 = []
lat2 = []
lat3 = []
labels = []
coordinates = []
limit = 0
projections = []
l_proj = []

## Read POSCAR. DO NOT MODIFY 
with open(INPUT_dir, 'r') as rf:
	for i, line in enumerate(rf):
		if i == 1:
			factor = float(line)
		if i == 2:
			for j in range(3):
				var =float(line.split()[j])
				lat1.append(factor*var)			
		if i == 3:
                       	for j in range(3):
                               	var =float(line.split()[j])
                               	lat2.append(factor*var)
		if i == 4:
                       	for j in range(3):
                               	var =float(line.split()[j])
                               	lat3.append(factor*var)
		if i == 5:
			lab = line.split()
		if i == 6:
			multiplicity = line.split()
			for j in range(len(multiplicity)):
				limit = limit + int(multiplicity[j])
		if i == 7:
			if line[0] == 'C' or line[0] == 'c':
				form = 'cartesian'
			if line[0] == 'D' or line[0] == 'd':
				form = 'direct'
		if i > 7 and i <= limit + 7:
			line = line.split()
			line_data = []
			for j in range(3):
				line_data.append(float(line[j]))
			coordinates.append(line_data)
my_el_pos = []

coun = 1
for i in range(len(multiplicity)):
	for j in range(int(multiplicity[i])):
		my_el_pos.append([coun,lab[i]])
		coun = coun + 1

A_mat = np.array([[lat1[0],lat2[0],lat3[0]],[lat1[1],lat2[1],lat3[1]],[lat1[2],lat2[2],lat3[2]]])
A_mat_inv = inv(A_mat)

res = []

if form == 'cartesian':
	for i in range(len(coordinates)):
		R_vec = np.array([float(coordinates[i][0]),float(coordinates[i][1]),float(coordinates[i][2])])
		res.append(np.dot(A_mat_inv,R_vec))

	for i in range(len(coordinates)):
		for j in range(len(coordinates[i])):
			coordinates[i][j] = res[i][j]

len_a, len_b, len_c = 0,0,0

for i in range(len(lat1)):
	len_a = len_a + lat1[i] * lat1[i]
len_a = math.sqrt(len_a)

for i in range(len(lat2)):
	len_b = len_b + lat2[i] * lat2[i]
len_b = math.sqrt(len_b)

for i in range(len(lat3)):
	len_c = len_c + lat3[i] * lat3[i]
len_c = math.sqrt(len_c)

cond_grid_x, cond_grid_y, cond_grid_z = 0,0,0

cond_grid_x = int(density/len_a)
cond_grid_y = int(density/len_b)
cond_grid_z = int(density/len_c)

if cond_grid_x%2==0:
	cond_grid_x = cond_grid_x + 1

if cond_grid_y%2==0:
        cond_grid_y = cond_grid_y + 1

if cond_grid_z%2==0:
        cond_grid_z = cond_grid_z + 1

with open(wannier90win_dir,'r+') as rf:
	for i, line in enumerate(rf):
		if line.strip() == "Begin Projections":
			start = i
		if line.strip() == "End Projections":
                        end = i

num = 0

with open(OUTCAR_dir,"r") as rf:
        for i, line in enumerate(rf):
                line = line.split()
                if line != [] and line[0] == "E-fermi":
                        efermi = float(line[2])

counter = 0
l_values = [[] for i in range(limit)]
my_el_w90 = []

with open(wannier90win_dir,'r+') as rf:
	for i, line in enumerate(rf):
		if i > start and i < end:
			counter = counter + 1
			line = line.replace(' = ','=')
			line = line.replace(':',' ')
			line = line.replace(';',' ')
			line = line.replace('!',' ')
			line = line.split()
			temp = []
			for i in range(1,len(line)):
				temp.append(line[i])
			my_el_w90.append(temp)

pos_list = []
for i in range(len(my_el_pos)):
	pos = -1
	for j in range(len(my_el_w90)):
		if int(my_el_pos[i][0]) == int(my_el_w90[j][len(my_el_w90[j])-2]) and str(my_el_pos[i][1]) == str(my_el_w90[j][len(my_el_w90[j])-1]):
			pos = j
	pos_list.append(pos)

for i in range(len(pos_list)):
	if pos_list[i] == -1:
		l_values[i].append(" ")
		projections.append(0)
	else:
		if "l=0" in my_el_w90[pos_list[i]]:
			l_values[i].append("l=0")
			num = num + 1
		if "l=1" in my_el_w90[pos_list[i]]:
			l_values[i].append("l=1")
			num = num + 3
		if "l=2" in my_el_w90[pos_list[i]]:
			l_values[i].append("l=2")
			num = num + 5

		projections.append(num)
		num = 0

projectors = [[] for i in range(limit)]

for j in range(limit):
	for i in l_values[j]:
		if i == 'l=0':
			projectors[j].append('s')
		if i == 'l=1':
			projectors[j].append('px')
			projectors[j].append('py')
			projectors[j].append('pz')
		if i == 'l=2':
			projectors[j].append('dxy')
			projectors[j].append('dyz')
			projectors[j].append('dzx')
			projectors[j].append('dx2-y2')
			projectors[j].append('dz2')
		if i == ' ':
			projectors[j].append(' ')

with open(wannier90wout_dir,'r') as rf:
	for i, line in enumerate(rf):
		line = line.split()
		if line and line[0] == "Final" and line[1] == "State":
			number = i
			break

wc = []
kt = []

if os.path.exists(kpts_band_dir):
	with open(kpts_band_dir,'r') as rf:
		for i, line in enumerate(rf):
			line = line.split()
			if line and i > 3:
				if line[4][0] == """\\""":
					kt.append([line[4][1],float(line[0]),float(line[1]),float(line[2])])
				else:
					kt.append([line[4],float(line[0]),float(line[1]),float(line[2])])

## Read wannier90.wout to extract wannier centers. DO NOT MODIFY
with open(wannier90wout_dir,'r') as rf:
	for i, line in enumerate(rf):
		line = line.replace(',',' ')
		line = line.replace('(','( ')
		line = line.replace(')',' )')
		line = line.split()
		if i > number and line and line[0] == "Sum":
                        break
		if i > number:
			line_data = []
			line_data.append(line[6])
			line_data.append(line[7])
			line_data.append(line[8])
			wc.append(line_data)

## Creation of submission file. MODIFY according to your needs
with open('run-wanntools.sh','w') as wf:
          wf.write("{}\n".format("#!/bin/bash"))
          wf.write("{}\n".format(f"#SBATCH -J {job_name}"))
          wf.write("{}\n".format(f"#SBATCH -A {project}"))
          wf.write("{}\n".format(f"#SBATCH -n {cores}"))
          wf.write("{}\n".format(f"#SBATCH --time={time}"))
          wf.write("{}\n".format("#SBATCH --export=ALL"))
          wf.write("{}\n".format(f"#SBATCH --mem-per-cpu={memory}"))
          wf.write("{}\n".format(f"#SBATCH -C {architecture}"))
          wf.write("\n")
          wf.write("{} {} {}\n".format("srun -n",cores,wanniertools_dir))

# Creation of wt.in file. MODIFY according to your needs
with open(OUTPUT_file,'w') as wf:
	wf.write("{}\n".format("&TB_FILE"))
	wf.write("{}\n".format("Hrfile = 'wannier90_hr.dat'"))
	wf.write("{}\n".format("Package = 'VASP'"))
	wf.write("{}\n".format("/"))
	wf.write("\n")
	wf.write("{}\n".format("LATTICE"))
	wf.write("{}\n".format("Angstrom"))
	for i in range(int(len(multiplicity))):
		for j in range(int(multiplicity[i])):
			labels.append(lab[i])

	for i in range(len(lat1)):
		wf.write("{} ".format(float(lat1[i])))
	wf.write("\n")
	for i in range(len(lat2)):
		wf.write("{} ".format(float(lat2[i])))
	wf.write("\n")
	for i in range(len(lat3)):
       		wf.write("{} ".format(float(lat3[i])))
	wf.write("\n")
	wf.write("\n")
	wf.write("{}\n".format("ATOM_POSITIONS"))
	wf.write("{}\n".format(limit))
	wf.write("{}\n".format("Direct"))
	for i in range(len(coordinates)):
		wf.write("{} {} {} {} \n".format(labels[i],coordinates[i][0],coordinates[i][1],coordinates[i][2]))

	wf.write("\n")
	wf.write("{}\n".format("PROJECTORS"))
	
	for i in range(limit):
		wf.write("{} ".format(projections[i]))
	wf.write("\n")
	for i in range(limit):
		wf.write("{} ".format(labels[i]))
		for j in projectors[i]:
			wf.write("{} ".format(j))
		wf.write("\n")
	wf.write("\n")
	wf.write("{}\n".format("SURFACE"))
	wf.write("{}\n".format("1  0  0"))
	wf.write("{}\n".format("0  1  0"))
	wf.write("\n")
	wf.write("{}\n".format("&CONTROL"))
	wf.write("{}\n".format("AHC_calc = T"))
	wf.write("{}\n".format("FindNodes_calc = F"))
	wf.write("{}\n".format("WeylChirality_calc = F"))
	wf.write("{}\n".format("BerryCurvature_calc = F"))
	wf.write("{}\n".format("BulkBand_calc = F"))
	wf.write("{}\n".format("BulkGap_Plane_calc = F"))
	wf.write("{}\n".format("BulkBand_calc = F"))
	wf.write("{}\n".format("BulkFS_calc = F"))
	wf.write("{}\n".format("BulkGap_cube_calc = F"))
	wf.write("{}\n".format("BulkGap_plane_calc = F"))
	wf.write("{}\n".format("SlabBand_calc = F"))
	wf.write("{}\n".format("WireBand_calc = F"))
	wf.write("{}\n".format("SlabSS_calc = F"))
	wf.write("{}\n".format("SlabArc_calc = F"))
	wf.write("{}\n".format("SlabQPI_calc = F"))
	wf.write("{}\n".format("SlabSpintexture_calc = F"))
	wf.write("{}\n".format("Wanniercenter_calc = F"))
	wf.write("{}\n".format("BerryCurvature_calc = F"))
	wf.write("{}\n".format("EffectiveMass_calc = F"))    
	wf.write("{}\n".format("/"))
	wf.write("\n")
	wf.write("{}\n".format("&SYSTEM"))
	wf.write("{}\n".format("SOC = 1"))
	wf.write("{}\n".format("NumOccupied = 2"))
	wf.write("{} {}\n".format("E_FERMI = ",efermi))
	wf.write("{}\n".format("/"))
	wf.write("\n")
	wf.write("{}\n".format("&PARAMETERS"))
	wf.write("{}\n".format("OmegaNum = 1001"))
	wf.write("{}\n".format("OmegaMin = -0.5"))
	wf.write("{}\n".format("OmegaMax = 0.5"))
	wf.write("Nk1 = {}\n".format(cond_grid_x))
	wf.write("Nk2 = {}\n".format(cond_grid_y))
	wf.write("Nk3 = {}\n".format(cond_grid_z))
	wf.write("{}\n".format("Gap_threshold = 0.0001"))
	wf.write("{}\n".format("/"))
	wf.write("\n")
	wf.write("{}\n".format("KCUBE_BULK"))
	wf.write("{}\n".format("0.00 0.00 0.00"))
	wf.write("{}\n".format("1.00 0.00 0.00"))
	wf.write("{}\n".format("0.00 1.00 0.00"))
	wf.write("{}\n".format("0.00 0.00 1.00"))
	wf.write("\n")
	wf.write("KPATH_BULK\n")
	wf.write("{}\n".format(int(len(kt)/2)))
	for i in range(len(kt)):
		wf.write("{} {} {} {} ".format(kt[i][0],kt[i][1],kt[i][2],kt[i][3]))
		if (i+1)%2 == 0:
			wf.write("\n")
	wf.write("\n")
	wf.write("{}\n".format("KPLANE_BULK"))
	wf.write("{}\n".format("0.00  0.00  0.00"))
	wf.write("{}\n".format("1.00  0.00  0.00"))
	wf.write("{}\n".format("0.00  1.00  0.00"))
	wf.write("\n")
    
    ## Wannier centers. DO NOT MODIFY
	wf.write("{}\n".format("WANNIER_CENTRES     ! copy from wannier90.wout"))
	wf.write("{}\n".format("Cartesian"))
	for i in range(len(wc)):
		wf.write("{} {} {}\n".format(wc[i][0],wc[i][1],wc[i][2]))

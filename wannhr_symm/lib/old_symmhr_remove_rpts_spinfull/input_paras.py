#!/usr/bin/env python
'''
set info of your crystal and wannier projections
'''
import sys, io
import numpy as np

# ------------------------------Change it according your case------------------------

# number of atoms in your crystal, consistent with your PROCAR
natom_all=4

# number of atoms taken into account in your wannier projections
natom=4

nptrans=1

# atom_num: No. of wannier's atoms which is consistent with that of OUTCAR's atoms' position
#                         part of CeCoIn5-OUTCAR 
#ion  position               nearest neighbor table
#  1  0.000  0.000  0.000-   4 3.25   4 3.25   4 3.25   4 3.25   6 3.28   3 3.28   5 3.28   6 3.28
#                            2 3.28   3 3.28   2 3.28   5 3.28
#  2  0.500  0.000  0.310-   7 2.71   7 2.71   3 2.86   6 3.25   6 3.25   6 3.25   6 3.25   1 3.28
#                            4 3.28   4 3.28   1 3.28
#  3  0.500  0.000  0.690-   7 2.71   7 2.71   2 2.86   5 3.25   5 3.25   5 3.25   5 3.25   1 3.28
#                            4 3.28   4 3.28   1 3.28
#  4  0.500  0.500  0.000-   1 3.25   1 3.25   1 3.25   1 3.25   6 3.28   2 3.28   2 3.28   3 3.28
#                            3 3.28   5 3.28   6 3.28   5 3.28
#  5  0.000  0.500  0.690-   7 2.71   7 2.71   6 2.86   3 3.25   3 3.25   3 3.25   3 3.25   4 3.28
#                            1 3.28   4 3.28   1 3.28
#  6  0.000  0.500  0.310-   7 2.71   7 2.71   5 2.86   2 3.25   2 3.25   2 3.25   2 3.25   4 3.28
#                            1 3.28   4 3.28   1 3.28
#  7  0.000  0.000  0.500-   2 2.71   3 2.71   5 2.71   6 2.71   5 2.71   6 2.71   2 2.71   3 2.71
#  .
#  .
#  . 
# irot  :   1
# --------------------------------------------------------------------
# isymop:   1   0   0
#           0   1   0
#           0   0   1

# gtrans:     0.0000000     0.0000000     0.0000000

# ptrans:     0.0000000     0.0000000     0.0000000
# rotmap:
# (   1->   1)  (   2->   2)  (   3->   3)  (   4->   4)  (   5->   5) 
# (   6->   6)  (   7->   7) 
#  
#  
# irot  :  2
# .
# .
#                       end of part of CeCoIn5-OUTCAR
# 
#                            part of wannier90.win                                 !   atom_num   
# begin projections
# f= 0.00000, 0.00000, 0.00000:l=3:z=0.0,0.0,1.0:x=0.5,0.5,0.0                     ! OUTCAR atom 1
# f= 0.00000, 0.00000, 0.00000:l=2:z=0.0,0.0,1.0:x=0.5,0.5,0.0                     !              
# f= 0.00000, 0.00000, 0.50000:l=2:z=0.0,0.0,1.0:x=1.0,0.0,0.0                     ! OUTCAR atom 7
# f= 0.00000, 0.50000, 0.31030:l=0;l=1:z=0.0,0.0, 1.0:x= 0.5,0.5,0.0               ! OUTCAR atom 6
# f= 0.50000, 0.00000, 0.31030:l=0;l=1:z=0.0,0.0, 1.0:x=-0.5,0.5,0.0               ! OUTCAR atom 2
# f= 0.00000, 0.50000, 0.68970:l=0;l=1:z=0.0,0.0,-1.0:x=-0.5,0.5,0.0               ! OUTCAR atom 5
# f= 0.50000, 0.00000, 0.68970:l=0;l=1:z=0.0,0.0,-1.0:x= 0.5,0.5,0.0               ! OUTCAR atom 3
# f= 0.50000, 0.50000, 0.00000:l=0;l=1:z=0.0,0.0, 1.0:x= 0.5,0.5,0.0               ! OUTCAR atom 4
# end projections
#                        end of part of wannier90.win
# sometimes atom_num may be a subset of all atoms. e.g. only OUTCAR'S atom-1 and 7, atom_num=[1,7]
# PLEASE check OUTCAR "  ion  position " block and wannier90.win projections, atoms_frac block 
atom_num =[1,2,3,4]
                     
# number of wannier orbitals
nband=16
norbs=32
nwann=32

# atoms position is same to wannier90.win not POSCAR 
# remeber set positions same both in OUTCAR and wannier90.win
# number of orbitals of each atom consistent with atom_num
#          Cefd  Ird In1p   ......  In5p
type_orbs=[['d'],['d'],['p'],['p']]
num_orbs =[[5],[5],[3],[3]]
lamda_orbs =[[0.00]]*4

# number of space group symmetry.
nsymm=8

# --------------------------End-Change it according your case------------------------

# check consistency
if len(num_orbs) != natom: 
   print "number of items in num_orbs should be same to natom in wannier projection"
   sys.exit(0)

if len(type_orbs) != natom: 
   print "number of items in type_orbs should be same to natom in wannier projection"
   sys.exit(0)

if len(atom_num) != natom: 
   print "number of items in atom_num should be same to natom in wannier projection"
   sys.exit(0)

# number of orbitals for each atom 
nors = np.zeros((natom),dtype=np.int32)

# index offset of each atom
offs = np.zeros((natom),dtype=np.int32)

# get offset and num-orbs, and print out infos
for jatom in range(natom):
   nors[jatom] = sum(num_orbs[jatom][:])
   offs[jatom] = sum(nors[0:jatom])

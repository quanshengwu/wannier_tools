#!/usr/bin/env python
'''
Define primitive lattice vectors, atoms' fractional coordinates and local 
axises for each atom
'''
from numpy import matrix,array,transpose,dot,cross,zeros
import numpy as np
from numpy.linalg import inv, det, norm
from math import sqrt, cos, sin
from input_paras import *

# define types and size
Amat=np.zeros((3,3),dtype=np.float64)
atoms_frac=zeros((3,natom),np.float64)
axis = zeros((3,3,natom),np.float64)
atoms_crac=zeros((3,natom),np.float64)

# ------------------------------Change it according your case------------------------
# primitive cell basis 
Amat[:,0]= [-1.7349141,  1.7349141,  5.8687908] # primitive cell basis 1
Amat[:,1]= [ 1.7349141, -1.7349141,  5.8687908] # primitive cell basis 2
Amat[:,2]= [ 1.7349141,  1.7349141, -5.8687908] # primitive cell basis 2
for i in range(3):
    print "len",i,norm(Amat[:,i])

# fractional axis in terms of primitive cell basis, of the atoms of WS supercell[000], 
# and it can be set in file wannier90.win projection block. When not set in this block,
# atoms_frac is of 1st cell with f^mu_i>=0 fraction coordinates of atoms in R=(0,0,0)
# WS super cell of Hamr.
atoms_frac[:,0]=  [0.833302,   0.333302,   0.50000] #La 1
atoms_frac[:,1]=  [0.083302,   0.083302,   0.00000] #La 2
atoms_frac[:,2]=  [0.251064,   0.751064,   0.50000] #La 2
atoms_frac[:,3]=  [0.501064,   0.501064,   0.00000] #La 2

# local axices for atoms in wannier projection, only set z and x axis, and
# y axis will be got by vector cross.
# Ce
for i in range(natom):
    axis[:,2,i] = [0.0,0.0,1.0]  # z 
    axis[:,0,i] = [1.0,0.0,0.0]  # x

# --------------------------End-Change it according your case------------------------

for iatom in range(natom): 
    for j in range(3): 
        if atoms_frac[j,iatom] < 0.0:
           atoms_frac[j,iatom] = atoms_frac[j,iatom] + 1.0
    print iatom+1, atoms_frac[:,iatom].T
                                           
for i in range(natom):
    print np.dot(Amat,atoms_frac[:,i]).T

for i in range(natom):
    axis[:,1,i] = cross(axis[:,2,i], axis[:,0,i]) # y
    for j in range(3):
         axis[:,j,i] = axis[:,j,i] / norm(axis[:,j,i]) # normalization, NECESSARY

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
Amat[:,0]= [     -3.8889403,    2.5564011,    4.1324340] # primitive cell basis 1
Amat[:,1]= [      3.8889403,   -2.5564011,    4.1324340] # primitive cell basis 2
Amat[:,2]= [      3.8889403,    2.5564011,   -4.1324340] # primitive cell basis 2

# fractional axis in terms of primitive cell basis, of the atoms of WS supercell[000], 
# and it can be set in file wannier90.win projection block. When not set in this block,
# atoms_frac is of 1st cell with f^mu_i>=0 fraction coordinates of atoms in R=(0,0,0)
# WS super cell of Hamr.
atoms_frac[:,0]=  [ 0.87156,  0.37613,  0.83685] #Sn   1
atoms_frac[:,1]=  [ 0.53927,  0.03471,  0.16315] #Te   1
atoms_frac[:,2]=  [ 0.46073,  0.62387,  0.49543] #Te   1
atoms_frac[:,3]=  [ 0.12844,  0.96529,  0.50457] #Te   1
atoms_frac[:,4]=  [ 0.12844,  0.62387,  0.16315] #Te   1
atoms_frac[:,5]=  [ 0.46073,  0.96529,  0.83685] #Te   1
atoms_frac[:,6]=  [ 0.53927,  0.37613,  0.50457] #Te   1
atoms_frac[:,7]=  [ 0.87156,  0.03471,  0.49543] #Te   1
atoms_frac[:,8]=  [ 0.92680,  0.10766,  0.81914] #Te   1
atoms_frac[:,9]=  [ 0.28852,  0.10766,  0.18086] #Te   1
atoms_frac[:,10]=  [ 0.71148,  0.89234,  0.81914] #Te   1
atoms_frac[:,11]=  [ 0.07320,  0.89234,  0.18086] #Te   1
atoms_frac[:,12]=  [ 0.20515,  0.03216,  0.82701] #Te   1
atoms_frac[:,13]=  [ 0.20515,  0.37814,  0.17299] #Te   1
atoms_frac[:,14]=  [ 0.79485,  0.96784,  0.17299] #Te   1
atoms_frac[:,15]=  [ 0.79485,  0.62186,  0.82701] #Te   1


#atoms_crac[:,0]= [-1.4913627,    4.4333000,    1.2340500] 
#atoms_crac[:,1]= [ 1.4913627,    4.4333000,    1.2340500]
#atoms_crac[:,2]= [-2.9265607,    6.1128114,    0.0000000]
#atoms_crac[:,3]= [ 2.9265607,    2.7537886,    0.0000000]
#atoms_crac[:,4]= [ 2.9265607,    6.1128114,    0.0000000]
#atoms_crac[:,5]= [-2.9265607,    2.7537886,    0.0000000]
#atoms_crac[:,6]= [ 0.6670235,    4.4333000,    0.0000000]
#atoms_crac[:,7]= [-0.6670235,    4.4333000,    0.0000000]
#atoms_crac[:,8]= [-2.4672919,    5.5989919,    1.2340500]
#atoms_crac[:,9]= [ 2.4672919,    3.2676081,    1.2340500]
#atoms_crac[:,10]=[ 2.4672919,    5.5989919,    1.2340500]
#atoms_crac[:,11]=[-2.4672919,    3.2676081,    1.2340500]
#
#invAmat = inv(Amat)
#for i in range(natom):
#    atoms_frac[:,i] = np.dot(invAmat,atoms_crac[:,i])
#    for j in range(3):
#        if atoms_frac[j,i] > 0.9999: atoms_frac[j,i] = atoms_frac[j,i] - 1.0
#        if abs(atoms_frac[j,i]) < 1.0E-4: atoms_frac[j,i] = 0.0



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

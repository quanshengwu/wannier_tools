#!/usr/bin/env python
import sys
import os
import numpy as np
from datetime import datetime

"""
Obtain a tight binding Hamiltonian of Haldane model with Wannier90 format

How to run

python haldane_hr_gen.py

This will generate the tight binding hamiltonian Haldane_hr.dat

LATTICE
Angstrom
2.1377110  -1.2342080   0.0000000
0.0000000   2.4684160   0.0000000
0.0000000   0.0000000   10.000000

ATOM_POSITIONS
2                               ! number of atoms for projectors
Direct                          ! Direct or Cartisen coordinate
C 0.333333 0.666667 0.500000
C 0.666667 0.333333 0.500000


"""
# Define tight-binding parameters
# You can find phase diagram in PRL 61,2015 (1988)
# Chern = 0
m=0.2; phi= np.pi/2.0; t1=1.0; t2=0.0;
# Gapless phase
#m=0.2; phi= np.pi/2.0; t1=1.0; t2=m/3.0/np.sqrt(3);
# Chern = 1
#m=0.2; phi= np.pi/2.0; t1=1.0; t2=m/3.0/np.sqrt(3)*2.0;


# maximum dimension for hr matrix
ndim = 2 
nrpts = 7
num_patom=2

# hr matrix 
norbs = num_patom*1
hmnr= np.zeros((norbs,norbs,nrpts),dtype = np.complex128)

# WS points
irvec   = np.zeros((3,nrpts),dtype = np.int32)

# degeneracy
dege   = np.zeros((nrpts),dtype = np.int32)+1


# complex unit
zi=1j

ir= 0
irvec[0, ir]= 0
irvec[1, ir]= 0
hmnr[0, 0, ir]=  m
hmnr[1, 1, ir]= -m
hmnr[0, 1, ir]=  t1
hmnr[1, 0, ir]=  t1

# 1 0 
ir= ir+1
irvec[0, ir]= 1
irvec[1, ir]= 0
hmnr[0, 0, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2
hmnr[0, 1, ir]= t1

# 0 1 
ir= ir+1
irvec[0, ir]= 0
irvec[1, ir]= 1
hmnr[0, 0, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2
hmnr[1, 0, ir]= t1

# 1 1 
ir= ir+1
irvec[0, ir]= 1
irvec[1, ir]= 1
hmnr[0, 0, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2

#-1 0 
ir= ir+1
irvec[0, ir]=-1
irvec[1, ir]= 0
hmnr[0, 0, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2
hmnr[1, 0, ir]= t1

# 0-1 
ir= ir+1
irvec[0, ir]= 0
irvec[1, ir]=-1
hmnr[0, 0, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2
hmnr[0, 1, ir]= t1

#-1-1 
ir= ir+1
irvec[0, ir]=-1
irvec[1, ir]=-1
hmnr[0, 0, ir]= (np.cos(phi)-zi*np.sin(phi)) *t2
hmnr[1, 1, ir]= (np.cos(phi)+zi*np.sin(phi)) *t2


#print "dump hr.dat..."
with open('Haldane_hr.dat','w') as f:
    line="Haldane model with m="+str(m)+", phi="+str(phi/np.pi)+"pi, t1="+str(t1)+", t2="+str(t2)+"Ref:Physical Review Letters 61, 18(1988)"+'\n'
    f.write(line)
    nl = np.int32(np.ceil(nrpts/15.0))
    f.write(str(norbs)+'\n')
    f.write(str(nrpts)+'\n')
    for l in range(nl):
        line="    "+'    '.join([str(np.int32(i)) for i in dege[l*15:(l+1)*15]])
        f.write(line)
        f.write('\n')
    for irpt in range(nrpts):
        rx = irvec[0,irpt];ry = irvec[1,irpt];rz = irvec[2,irpt]
        for jatomorb in range(norbs):
            for iatomorb in range(norbs):
               rp =hmnr[iatomorb,jatomorb,irpt].real
               ip =hmnr[iatomorb,jatomorb,irpt].imag
               line="{:8d}{:8d}{:8d}{:8d}{:8d}{:20.10f}{:20.10f}\n".format(rx,ry,rz,jatomorb+1,iatomorb+1,rp,ip)	
               f.write(line)

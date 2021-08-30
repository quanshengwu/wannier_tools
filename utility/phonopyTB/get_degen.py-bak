#!/usr/bin/env python
import sys
import os
import numpy as np
#from input_vectors import *

# Amat(:,i): a_i, primitive cell vectors
def get_rpts(supercell_dim1,supercell_dim2,supercell_dim3):
    ndiff = np.zeros((3),dtype=np.int32)
    ndegen = np.zeros((10000),dtype=np.int32)
    dist = np.zeros((125),dtype=np.float64)
    irvec = np.zeros((10000,3),dtype=np.float64)
    
    mp_grid =[supercell_dim1,supercell_dim2,supercell_dim3] 
    nrpts = 0
    for n1 in range(-mp_grid[0],mp_grid[0]+1):
       for n2 in range(-mp_grid[1],mp_grid[1]+1):
           for n3 in range(-mp_grid[2],mp_grid[2]+1):
                 print "------------n1 n2 n3", n1, n2, n3
                 icnt = 0  
                 rpts=[]
                 for i1 in range(-2,2+1):
                    for i2 in range(-2,2+1):  
                       for i3 in range(-2,2+1):  
                          # Calculate distance squared |r-R|^2
                          ndiff[0] = n1 - i1 * mp_grid[0]  
                          ndiff[1] = n2 - i2 * mp_grid[1]  
                          ndiff[2] = n3 - i3 * mp_grid[2]  
                          dist[icnt] = np.linalg.norm(np.dot(Amat,ndiff))
                          rpts.append([icnt,ndiff[0],ndiff[1],ndiff[2]])
                          icnt = icnt + 1 
                 dist_min=min(dist)
                 if abs(dist[62] - dist_min) < 1.0E-7:
                    ndegen[nrpts]=0
                    for i in range(125):
                       if (abs (dist [i] - dist_min) < 1.0E-7 ): ndegen[nrpts]=ndegen[nrpts]+1
                       if (abs (dist [i] - dist_min) < 1.0E-7 ): print rpts[i]
                    irvec[nrpts,0] = n1  
                    irvec[nrpts,1] = n2   
                    irvec[nrpts,2] = n3
                    if (n1==0  and  n2==0  and  n3==0): rpt_origin=nrpts
                    print nrpts, "|", irvec[nrpts,:], "|",ndegen[nrpts] 
                    nrpts = nrpts+1
    # Check the "sum rule"
    tot = 0.0
    for i in range(nrpts):
       tot = tot + 1.0/np.real(ndegen[i]) 
    print "tot sum rule", tot
    if (abs (tot - np.real(mp_grid[0] * mp_grid[1] * mp_grid[2]) ) > 1.0E-8):
       print "sum rule wrong!"
    degen = ndegen[0:nrpts]
    rpts = irvec[0:nrpts]
    return nrpts, degen, rpts 
    
nrpt,rpts,ndegen=get_rpts(3,3,3)
print nrpt
print ndegen
# get 

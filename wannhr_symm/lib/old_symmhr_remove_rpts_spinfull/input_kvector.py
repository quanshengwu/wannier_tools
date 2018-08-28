#!/usr/bin/env python
'''
Developer: Changming Yue, yuechangming8@gmail.com
set info of your reciprocal basis and high symm k.
'''
import sys, io
import numpy as np
from input_vectors import Amat

# ------------------------------Change it according your case------------------------
# high symmetry K-point in 1st BZ.
hsymkpt=[[0.500,0.250,0.750],# W
         [0.500,0.500,0.500],# L
         [0.000,0.000,0.000],# G
         [0.500,0.000,0.500],# X
         [0.500,0.250,0.750],# W
         [0.375,0.375,0.750]]# K

# number of kpoints on a k-path, say from G-X
nkpt_path = 50
# --------------------------End-Change it according your case------------------------
 
# more precise than kbase by direct calc from lattice vectors
kbase=np.zeros((3,3),dtype=np.float64)
vol = np.dot(Amat[:,0],np.cross(Amat[:,1],Amat[:,2]))
kbase[:,0] = np.cross(Amat[:,1],Amat[:,2])/vol
kbase[:,1] = np.cross(Amat[:,2],Amat[:,0])/vol
kbase[:,2] = np.cross(Amat[:,0],Amat[:,1])/vol



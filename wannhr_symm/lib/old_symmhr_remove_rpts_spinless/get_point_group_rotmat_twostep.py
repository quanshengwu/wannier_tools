#!/usr/bin/env python
'''
get rotation matrices of all symmetry opreations defined in symmop.dat file.
'''
from numpy import matrix, complex,array,transpose,dot,cross,zeros,reshape,transpose
from numpy.linalg import inv, det,eigh, norm
import numpy as np
from math import sqrt, cos, sin
from input_vectors import *
from input_paras import *

from get_orb_rotmat_twostep import get_any_rot_orb_twostep
from read_symop import *
import wanntb.rotate as rot_spin
from get_euler_angle import *

# rot matrix in terms of primitive cell basis
rot = zeros((3,3),np.float64)

# rot matrix in terms of cartecian axis
rot_glb = zeros((3,3),np.float64)

# rot matrix in terms of local axis
rot_loc = zeros((3,3),np.float64)

gtrans = zeros((3),np.float64)
ptrans = zeros((3),np.float64)
new_vec = zeros((3),np.float64)
old_vec = zeros((3),np.float64)

# rot map for atoms
rotmap = zeros((natom),np.int32)

# all orbitals' rot map 
vec_shift  = zeros((nsymm,nptrans,natom,3),np.float64)

nors = np.zeros((natom),dtype=np.int32)
offs = np.zeros((natom),dtype=np.int32)
# get off set of each atom. get number of orbitals 
for jatom in range(natom):
   nors[jatom] = sum(num_orbs[jatom][:])
   offs[jatom] = sum(nors[0:jatom])
   print "atom", jatom+1, "off", offs[jatom]+1,  "norb", nors[jatom]
    
invAmat=inv(Amat)

# rot matrix: {\hat P}_{R} |j,atom B>=\sum_{i}P_{R|i,j}|i,atom B'>
prottmp=np.zeros((nband,nband,nptrans,nsymm),dtype=np.complex128)
protmat=np.zeros((norbs,norbs,nptrans,nsymm),dtype=np.complex128)
for isymm in range(0,nsymm):
    print "============================================symm i=",isymm+1
    for ii in range(0,3):
        for jj in range(0,3):
            rot[ii,jj] = np.float64(symop[isymm][0][ii][jj])

    print "rot\n", rot

# if avec is written in A=(a1,a2,a3) 
#  a1  a2  a3
# a11 a12 a13
# a21 a22 a23 
# a31 a32 a33
# then rot_glb=R'=A.R.A^-1. if you act R' on an vector written in cartetian axis, you
# get new vector also written in cartician axis
    cmat = dot(Amat,rot)
    rot_glb = dot(cmat,invAmat)

# rot_glb often is 1 or -1. if too small, it may arise from number error, and force it
# to be zero in this case
    for i in range(3): 
        for j in range(3):
            if abs(rot_glb[j,i])<1E-6:
                rot_glb[j,i]=0.0

    print "rot_glb=\n", rot_glb
    print "det of rot_glb=",det(rot_glb)

    for ii in range(0,3):
        gtrans[ii] = np.float64(symop[isymm][1][ii])
    print "gtrans\n", gtrans

    for iptr in range(nptrans):
        for ii in range(0,3):
             ptrans[ii] = np.float64(symop[isymm][3][iptr][ii])
        print "ptrans", iptr+1, "\n", ptrans

        for ii in range(0,natom):
            rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]
    
        print "atom rotmap=\n", 
        for i in range(natom):
            print i+1, "->", rotmap[i]+1
   
    # loop over all atoms, determine each one's orbital rot map
        offj = 0
        for jatom in range(0,natom):
            iatom = rotmap[jatom]
            print "--------------------------------"
            print "jatom->iatom", jatom+1, iatom+1
    
            old_vec = atoms_frac[:,jatom]
            new_vec = dot(rot,old_vec) + gtrans + ptrans
            vec_shift[isymm,iptr,jatom] = new_vec - atoms_frac[:,iatom]

            # round it to make integers
            for jj in range(3):
                vec_shift[isymm,iptr,jatom,jj] = np.round(vec_shift[isymm,iptr,jatom,jj])
            print "isymm", isymm+1, "iptr",iptr+1, "jatom", jatom+1, "vec_shift", vec_shift[isymm,iptr,jatom]
    
           #G_vec[isymm][jatom] = new_vec + gtrans - atoms_frac[:,iatom]
            axs_loc_jatom = axis[:,:,jatom]
            axs_loc_iatom = axis[:,:,iatom]
    
    # loop over jatom's all kind of orbitals. e.g. Ce-f,d
            for k in range(len(type_orbs[jatom])):
                norb = num_orbs[jatom][k] 
                print "jatom",jatom+1,"orbital",type_orbs[jatom][k], "offj", offj+1, "norb", norb, "-------"
    
    # first step, P_{rot_glb} |phi a jatom> = \sum_a' P_{a'a} | phi a' iatom> with jatom local axis
                rot_axs1 = dot(dot(inv(axs_loc_jatom),rot_glb),axs_loc_jatom)
                #rot_orb1 = get_any_rot_orb_twostep(type_orbs[jatom][k],rot_axs1)
                rot_axs2 = dot(transpose(axs_loc_iatom),axs_loc_jatom)
                #rot_orb2 = get_any_rot_orb_twostep(type_orbs[jatom][k],rot_axs2)
                rot_axs = dot(rot_axs2,rot_axs1)
                #rot_orb  = np.dot(rot_orb2,rot_orb1)
                rot_orb  = get_any_rot_orb_twostep(type_orbs[jatom][k],rot_axs)
    
                print "rot_orb\n", rot_orb 
                prottmp[offj:offj+norb,offj:offj+norb,iptr,isymm] = rot_orb
                offj = offj + norb
    # P_{rot} for spinors: up dn up dn ...
        protmat[:,:,iptr,isymm] = prottmp[:,:,iptr,isymm]
    #   for i in range(norbs):
    #       for j in range(norbs):
    #           if (abs(protmat[j,i,isymm])>1.0E-8): print "{:6d}{:6d}{:12.8f}{:12.8f}".format(j+1, i+1, protmat[j,i,isymm].real, protmat[j,i,isymm].imag)
        print "+++++++++++++++++++++++++++++++++++++++++++++++++"    

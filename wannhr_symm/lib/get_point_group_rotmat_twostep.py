#!/usr/bin/env python
'''
get rotation matrices of all symmetry opreations defined in symmop.dat file.
'''
from numpy import matrix, complex,array,transpose,dot,cross,zeros,reshape,transpose
from numpy.linalg import inv, det,eigh, norm
import numpy as np
from math import sqrt, cos, sin

from get_orb_rotmat_twostep import get_any_rot_orb_twostep
import rotate as rot_spin
from get_euler_angle import *
#from dets_rot import rot_glb_dets
import os

np.set_printoptions(precision=10,linewidth=200)

lst_prec = np.zeros((3),dtype=np.float64) 
hsq2 = 0.5 * 2.0**0.5
hsq3 = 0.5 * 3.0**0.5
hsq1 = 0.5
lst_prec[0]=hsq1
lst_prec[1]=hsq2
lst_prec[2]=hsq3

# get rotation matrix in terms of all projected wannier orbitals
def get_rotation_matrix_in_wannorbs(myposwan, nsymm, nptrans, symop, wann_atom_rotmap):
    Amat        = myposwan.Amat
    natom_pos   = myposwan.natom_pos     
    atoms_frac  = myposwan.atoms_frac    
    atoms_cart  = myposwan.atoms_cart    
    atoms_symbol= myposwan.atoms_symbol  
    ispinor     = myposwan.ispinor       
    natom_wan   = myposwan.natom_wan     
    type_orbs   = myposwan.type_orbs     
    num_orbs    = myposwan.num_orbs      
    atom_num    = myposwan.atom_num      
    axis        = myposwan.axis          
    nband       = myposwan.nband
    norbs       = myposwan.norbs

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
    rotmap = zeros((natom_wan),np.int32)
    
    # all orbitals' rot map 
    vec_shift  = zeros((nsymm,nptrans,natom_wan,3),np.float64)
    
    nors = np.zeros((natom_wan),dtype=np.int32)
    offs = np.zeros((natom_wan),dtype=np.int32)
    # get off set of each atom. get number of orbitals 
    print "----------To get rot mat(protmat)in wann orbs------"
    print "########See ./protmat.log file for details#########"
    print "" 
    f = open('protmat.log', 'w')
    print >>f, "offset of each atoms in protmat:"
    for jatom in range(natom_wan):
       nors[jatom] = sum(num_orbs[jatom][:])
       offs[jatom] = sum(nors[0:jatom])
       print >>f, "atom", jatom+1, "off", offs[jatom]+1,  "norb", nors[jatom]
        
    print >>f, ""
    invAmat=inv(Amat)
    
    # rot matrix: {\hat P}_{R} |j,atom B>=\sum_{i}P_{R|i,j}|i,atom B'>
    prottmp=np.zeros((nband,nband,nptrans,nsymm),dtype=np.complex128)
    protmat=np.zeros((norbs,norbs,nptrans,nsymm),dtype=np.complex128)
    dmat   =np.zeros((2,2,nsymm),dtype=np.complex128)
    for isymm in range(nsymm):
        if os.path.exists("write_symmop.flag"): f2 = open('Representation_Matrix_symop_'+str(isymm+1)+'.dat', 'w')
	if os.path.exists("write_symmop.flag"): f3 = open('dmat_symop_'+str(isymm+1)+'.dat', 'w')
        print >>f, "============================================symm i=",isymm+1
        rot = symop[isymm][0]
        print >>f, "rot(Bravis)\n", rot
    
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
        print >>f, "rot_glb before forecing symm =\n", rot_glb	
	rgb = rot_glb * 1.0
	for i in range(3):
	    for j in range(3):
	        a = rot_glb[j,i]
		if abs(a) < 0.01: 
		    rot_glb[j,i] = 0.0
		else: 
		    for k in range(3):
		        m = lst_prec[k]
                        if abs(abs(a)-m)<0.01 : rot_glb[j,i] = (a/abs(a)) * m


        print >>f, "rot_glb after  forecing symm =\n", rot_glb
	print >>f, "diffs\n", rgb - rot_glb

        print >>f, "det of rot_glb=",det(rot_glb)
        #rot_glb_dets.append(det(rot_glb))
    
    # get dmat
        if ispinor: 
            euler_ang = np.zeros((3),dtype=np.float64)
            euler_ang[0],euler_ang[1],euler_ang[2] = rmat2euler(det(rot_glb)*rot_glb)
            dmat[:,:,isymm] = rot_spin.dmat_spinor(euler_ang[0],euler_ang[1],euler_ang[2])
	    print>>f, "euler angles pi", euler_ang/np.pi
            tmp = rot_spin.euler_to_rmat(euler_ang[0],euler_ang[1],euler_ang[2])
	    print>>f, "error of rot_glb and that by euler2rmat\n", rot_glb - det(rot_glb)*tmp
	    print>>f, "euler2rmat sum abs error:", sum(sum(abs(rot_glb-det(rot_glb)*tmp)))

	    #dmat[:,:,isymm]= spin_reps(rot_glb) # G Winker's 
         
            print>>f, "dmat before forcing symm = \n", dmat[:,:,isymm]
            # correct numerical errors
	    dtmp = dmat[:,:,isymm] * 1.0
            for i in range(2):
                for j in range(2):
		    a = dmat[j,i,isymm]
		    if abs(a) < 0.01: 
		        dmat[j,i,isymm] = 0.0
		    else: 
		        for k in range(3):
		            m = lst_prec[k]
                            if abs(abs(a)-m)<0.0001 : dmat[j,i,isymm] = (a/abs(a)) * m

            print>>f, "dmat after  forcing symm = \n", dmat[:,:,isymm]
	    print>>f, "diffs\n", dtmp - dmat[:,:,isymm]
     
        gtrans = np.float64(symop[isymm][1])
        print >>f, "gtrans\n", gtrans
    
        for iptr in range(nptrans):
            ptrans = np.float64(symop[isymm][3][iptr])
            print >>f, "ptrans", iptr+1, "\n", ptrans
    
            for ii in range(0,natom_wan):
                rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]
        
            print >>f, "wann atom rotmap=\n", 
            for i in range(natom_wan):
                print >>f, i+1, "->", rotmap[i]+1
       
        # loop over all atoms, determine each one's orbital rot map
            offj = 0
            for jatom in range(0,natom_wan):
                iatom = rotmap[jatom]
                print >>f, "--------------------------------"
                print >>f, "jatom->iatom", jatom+1, iatom+1
        
                old_vec = atoms_frac[:,atom_num[jatom]-1]
                new_vec = dot(rot,old_vec) + gtrans + ptrans
                print >>f, "old_vec = atoms_frac[:,atom_num[jatom]-1]   ", old_vec
                print >>f, "new_vec = dot(rot,old_vec) + gtrans + ptrans", new_vec
                print >>f, "          atoms_frac[:,atom_num[iatom]-1]   ", atoms_frac[:,atom_num[iatom]-1]
                vec_shift[isymm,iptr,jatom] = new_vec - atoms_frac[:,atom_num[iatom]-1]
    
                # round it to make integers
                for jj in range(3):
                    vec_shift[isymm,iptr,jatom,jj] = np.round(vec_shift[isymm,iptr,jatom,jj])
                print >>f, "isymm", isymm+1, "iptr",iptr+1, "jatom", jatom+1, "vec_shift            ", vec_shift[isymm,iptr,jatom]
        
               #G_vec[isymm][jatom] = new_vec + gtrans - atoms_frac[:,iatom]
                axs_loc_jatom = axis[:,:,jatom]
                axs_loc_iatom = axis[:,:,iatom]
        
        # loop over jatom's all kind of orbitals. e.g. Ce-f,d
                for k in range(len(type_orbs[jatom])):
                    norb = num_orbs[jatom][k] 
                    print >>f, "jatom",jatom+1,"orbital",type_orbs[jatom][k], "offj", offj+1, "norb", norb, "-------"
        
        # first step, P_{rot_glb} |phi a jatom> = \sum_a' P_{a'a} | phi a' iatom> with jatom local axis
                    rot_axs1 = dot(dot(inv(axs_loc_jatom),rot_glb),axs_loc_jatom)
                    rot_axs2 = dot(transpose(axs_loc_iatom),axs_loc_jatom)
                    rot_axs = dot(rot_axs2,rot_axs1)
                    rot_orb  = get_any_rot_orb_twostep(type_orbs[jatom][k],rot_axs)
                    #print  "FUCK ALL T2G", type_orbs[jatom][k]
                    #rot_orb  = get_any_rot_orb_twostep('t2g',rot_axs)
        
                    print >>f, "protrotmat_orb\n", rot_orb 
                    prottmp[offj:offj+norb,offj:offj+norb,iptr,isymm] = rot_orb
                    offj = offj + norb
            if ispinor:
                # P_{rot} for spinors: up dn up dn ...
                protmat[:,:,iptr,isymm] = np.kron(prottmp[:,:,iptr,isymm],dmat[:,:,isymm])
            else:
                protmat[:,:,iptr,isymm] = prottmp[:,:,iptr,isymm]
            #for i in range(norbs):
            #    for j in range(norbs):
            #        if (abs(protmat[j,i,isymm])>1.0E-8): print "{:6d}{:6d}{:12.8f}{:12.8f}".format(j+1, i+1, protmat[j,i,isymm].real, protmat[j,i,isymm].imag)
            print >>f, "+++++++++++++++++++++++++++++++++++++++++++++++++"    
        for i1 in range(2):
            for i2 in range(2):
	        line="{:8d}{:8d}{:24.14f}{:24.14f}".format(i1,i2,dmat[i1,i2,isymm].real, dmat[i1,i2,isymm].imag) 
		if os.path.exists("write_symmop.flag"): print >>f3, line
        for i1 in range(nband):
            for i2 in range(nband):
		if abs(prottmp[i1,i2,iptr,isymm])>1.0E-6:
 	            line="{:8d}{:8d}{:24.14f}{:24.14f}".format(i1,i2,prottmp[i1,i2,iptr,isymm].real, prottmp[i1,i2,iptr,isymm].imag)
  		    if os.path.exists("write_symmop.flag"): print >>f2, line
        if os.path.exists("write_symmop.flag"): f2.close()
        if os.path.exists("write_symmop.flag"): f3.close()
        
    f.close() 
    return nors, offs, vec_shift, protmat

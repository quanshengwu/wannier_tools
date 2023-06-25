#!/usr/bin/evn python

import sys
import spglib
import numpy as np
from atoms import Atoms

def show_symmetry(symmetry):
    for i in range(symmetry['rotations'].shape[0]):
        print("----------------- %4d ----------------------------" % (i + 1))
        rot = symmetry['rotations'][i]
        trans = symmetry['translations'][i]
        print("  rotation:")
        for x in rot:
            print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
        print("  translation:")
        print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))

def show_lattice(lattice):
    print "-------------------poscar--------------------------"
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        lth = np.dot(vec,vec)**0.5
        print("%s %12.8f %12.8f %12.8f" % (tuple(axis,) + tuple(vec)))

def show_cell(lattice, positions,labels):
    show_lattice(lattice)
    print("Atomic points:")
    i = 1
    for p, s in zip(positions, labels):
        print("%4s %4s %10.5f %10.5f %10.5f" % ((str(i),)+(s,)+tuple(p)))
        i = i+1

def get_rot_trans(myposwan):
    Amat        = myposwan.Amat
    natom_pos   = myposwan.natom_pos     
    atoms_frac  = myposwan.atoms_frac    
    atoms_symbol= myposwan.atoms_symbol  
    natom_wan   = myposwan.natom_wan     
    atom_num    = myposwan.atom_num      

    crystal_ase = Atoms(symbols=atoms_symbol,cell=Amat.T, scaled_positions=atoms_frac.T,pbc=True)
   
    symmetry = spglib.get_symmetry(crystal_ase)
    ptgrp = spglib.get_pointgroup(symmetry['rotations'])
    
    nsymm = symmetry['rotations'].shape[0]
    
    rotmat = symmetry['rotations']
    rottrn = symmetry['translations']
    rmat0 = np.array(rotmat[0], dtype = np.float64)
    
    # get the number of ptrans
    nptrans = 0
    for i in range(nsymm):
        rmati = np.array(rotmat[i], dtype = np.float64)
        diff = abs(rmati-rmat0).sum()
        #print "i ", i+1, "diff", diff
        if diff < 1.0E-9: nptrans = nptrans + 1
    
    # get the number of distinct rotaions
    nrot = nsymm/nptrans

    print "--------------abstrac of symmetry------------------"
    print("Spacegroup of this cell is %s." % spglib.get_spacegroup(crystal_ase))
    print('Point group: ',ptgrp[0])
    print "number of total symmetry operations: ", nsymm
    print "number of distnct rot  ", nrot 
    print "nptrans per rotation", nptrans
    print "######### For details, see ./mysymmop.dat ########"
    print ""
    
    # check whether rotmat are put as typeI
    # still need to be tested. 
    err="""
    Opps, It seems that in spglib, the rotaions_matrix is not put as
    ----typeI-----
    rot1/ptran_11
    rot2/ptran_21  
    ...
    rotn/ptran_n1  
    ...
    ...
    ...
    rot1/ptran_1m 
    rot2/ptran_2m
    ...
    rotn/ptran_nm
    
    If put as above, I want to trans it as in vasp/OUTCAR
    
    ----typeII-----
    rot1/ptran_11
    rot1/ptran_12
    ...
    rot1/ptran_1m
    ...
    ...
    ...
    rotn/ptran_n1  
    rotn/ptran_n2  
    ...
    rotn/ptran_nm  
    
    """
    for isym in range(nrot):
        rmat0 = np.array(rotmat[isym], dtype = np.float64)
        diff = 0.0
        for iptr in range(nptrans): 
            numb = isym + iptr*nrot
            rmatt = np.array(rotmat[numb], dtype = np.float64)
            diff = diff + abs(rmatt-rmat0).sum() 
        if diff>1.0E-5: print err
        if diff>1.0E-5: sys.exit(0)
    
    # To change code less, we use type I and set nptrans=1
    nptrans = 1 
    wann_atom_rotmap=np.zeros((nsymm,nptrans,natom_wan), dtype=np.int32)
    symop = []
    f = open('mysymmop.dat', 'w')
    for isym in range(nsymm):
        print>>f, ("  --------------- %4d ---------------" % (isym + 1))
        rot = symmetry['rotations'][isym]
        trans = symmetry['translations'][isym]
        print>>f, ("  rotation:")
        for x in rot:
            print>>f, ("      %2d %2d %2d " % (x[0], x[1], x[2]))
        print>>f, ("  translation:")
        print>>f, ("      %8.5f %8.5f %8.5f " % (trans[0], trans[1], trans[2]))
        print>>f, ("  rotmap:")
        gtrans = trans
        ptrans = [np.zeros((3),dtype=np.float64)]
# rot map for POSCAR
        rotmap = []
        for iatom in range(natom_pos):    
            atoms_frac_new = np.dot(rotmat[isym],atoms_frac[:,iatom]) + rottrn[isym]
            # find equiv loc as in original cell
            for jatom in range(natom_pos):
               frac_diff = atoms_frac[:,jatom]-atoms_frac_new
	       # determine the equivalent between two atoms before and after space group operation.
               if abs(frac_diff-np.array([round(i) for i in frac_diff],dtype=np.float64)).sum() < 1.0E-2: 
                   print>>f, "             ", iatom+1, "->", jatom+1
                   rotmap.append([iatom+1,jatom+1])
	       
# rot map for wann.in. here 
        for iatom in range(natom_wan):    
            target = rotmap[ atom_num[iatom]-1 ][1]
            wann_atom_rotmap[isym,0,iatom]=atom_num.index(target)
        symop.append((rotmat[isym],gtrans,rotmap,ptrans))
        print>>f, "\n"
        print>>f, "\n"

    print>>f,"rotations as for atoms in wannier projection" 
    # rot map of atoms in wannier projection
    for isymm in range(nsymm):
        print >>f, "isymm", isymm+1, "---------"
        for iptr in range(nptrans): 
            print>>f, "trans", 0, gtrans
            for ii in range(natom_wan): 
                print>>f, "OUTCAR atom rot map:", atom_num[ii], "->", symop[isymm][2][atom_num[ii]-1][1], "=====>", "wannier atom rot map:", ii+1, "->", wann_atom_rotmap[isymm,iptr,ii]+1
    f.close()
    return nsymm, nptrans, symop, wann_atom_rotmap

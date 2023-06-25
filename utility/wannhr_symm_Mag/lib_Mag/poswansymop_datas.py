#!/usr/bin/env python
import sys

import numpy as np
from read_poswan import poswan
from get_symmop import get_rot_trans, show_cell
from get_point_group_rotmat_twostep import get_rotation_matrix_in_wannorbs
myposwan    = poswan.get_poswan()
Amat        = myposwan.Amat
natom_pos   = myposwan.natom_pos     
atoms_frac  = myposwan.atoms_frac    
atoms_cart  = myposwan.atoms_cart    
atoms_magmom= myposwan.atoms_magmom
atoms_symbol= myposwan.atoms_symbol  
ispinor     = myposwan.ispinor       
natom_wan   = myposwan.natom_wan     
type_orbs   = myposwan.type_orbs     
num_orbs    = myposwan.num_orbs      
atom_num    = myposwan.atom_num      
axis        = myposwan.axis          
nband       = myposwan.nband
norbs       = myposwan.norbs

# print out pocar
show_cell(Amat.T,atoms_frac.T,atoms_symbol)

print "---------------------------------------------------"
print "Length of Bravis vectors. Maybe you need to check it"
print "a1:", np.dot(Amat[:,0],Amat[:,0])**0.5
print "a2:", np.dot(Amat[:,1],Amat[:,1])**0.5
print "a3:", np.dot(Amat[:,2],Amat[:,2])**0.5
print ""                                           
 
print "-------------------atoms cart----------------------"
for iatom in range(natom_pos):
    a, b, c = atoms_cart[:,iatom]
    print("%4d %2s %12.8f %12.8f %12.8f" % (iatom+1, atoms_symbol[iatom], a,b,c))
print ""                                           

print "----------------orbital projected------------------"
print "Make sure that the centers of wannier orbitals are consistent with those shown in wannier90.wout/spread"
for iatom in range(natom_wan):
    a, b, c = atoms_cart[:,atom_num[iatom]-1]
    print "wann No.", iatom+1, "poscar No.", atom_num[iatom], " ", atoms_symbol[atom_num[iatom]-1], " orbs", type_orbs[iatom], "dim",  num_orbs[iatom], "cart", ("{:12.8f}"*3).format(a,b,c)
print ""
print "spinor? ", ispinor
print ""
print "---------------number of oritals-------------------" 
print "nband=", nband, "norbs=", norbs 
print ""

nsymm, nptrans, symop, wann_atom_rotmap = get_rot_trans(myposwan)
nors, offs, vec_shift, protmat = get_rotation_matrix_in_wannorbs(myposwan, nsymm, nptrans, symop, wann_atom_rotmap)

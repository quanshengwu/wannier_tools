#!/usr/bin/env python
'''
Developer: Changming Yue, yuechangming10@gmail.com
perform symmetrization for tight-binding Hamiltonian abtained by 
wannier110. For equivalent atoms, local axices should be set in
wannier110.win file. 
'''
import datetime,sys,io
import numpy as np

from poswansymop_datas import *
from read_hamr import HR

np.set_printoptions(precision=6,suppress=True,linewidth=150)
sdim = 12
# rot matrix in terms of primitive cell basis
rot = np.zeros((3,3),np.float64)
rotmap = np.zeros((natom_wan),np.int32) # rot map for atoms
rpt = np.zeros((3), dtype=np.float64)

with open('rpt_hop.dat', 'r') as f:
     a,b,c,d,e = f.readline().strip().split()
     rpt[0] = np.float64(a)
     rpt[1] = np.float64(b)
     rpt[2] = np.float64(c)
     iatom = np.int32(d)
     jatom = np.int32(e)
     print a, b, c, "|", d, e

offj=2*offs[jatom] # off set of jatom in norbs 
norj=2*nors[jatom] # number of orbitals of jatom
offi=2*offs[iatom]
nori=2*nors[iatom]


f = open('rpt_map.log', 'w')
for isymm in range(nsymm):         
    rot = np.float64(symop[isymm][0])
    rot_rpt = np.dot(rot,rpt) 
    for iptr in range(nptrans):
        for ii in range(0,natom_wan):
            rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]
        offjp=2*offs[rotmap[jatom]] # off set of the atom which jatom rot to in norbs
        norjp=2*nors[rotmap[jatom]] # number of orbitals of the atom which jatom rot to 
        offip=2*offs[rotmap[iatom]]
        norip=2*nors[rotmap[iatom]]


        new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
        ixt = np.int32(np.round(new_rpt[0])) 
        iyt = np.int32(np.round(new_rpt[1])) 
        izt = np.int32(np.round(new_rpt[2]))
        if abs(ixt)>sdim or abs(iyt)>sdim or abs(izt)>sdim: # for my case, ixyzt>10 may occur once in nsymm. so here should dnorm-1
            print "--------------------Out of Huge Box, irtp", rpt, "-->", new_rpt, "-------------------"
	    print "Error! please open symmhrspinless.py and make a bigger sdim!!!!!! "
	    sys.exit(0)
        else: 
	    # data at new rpts should be otained by known datas  
            print>>f, "            hr.rpts:", rpt, "|", iatom, jatom, "-->", new_rpt, "|", rotmap[iatom], rotmap[jatom ], "|", offi+1, offj+1, "->", offip+1, offjp+1


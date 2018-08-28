#!/usr/bin/env python
'''
Developer: Changming Yue, yuechangming10@gmail.com
perform symmetrization for tight-binding Hamiltonian abtained by 
wanniersdim+sdim. For equivalent atoms, local axices should be set in
wanniersdim+sdim.win file. 
'''
import datetime,sys,io
import numpy as np

from input_vectors import *
from input_paras import *
from input_kvector import *
from wanntb.kvec import SymKVec, UniKVec
from wanntb.tran import fourier_hr2hk
from get_point_group_rotmat_twostep import *
from read_symop import *

from numpy.linalg import eigh 
from wanntb.ham import HR
from wanntb.io import write_band
from progressbar import *

def get_equiv_rpts(mp_grid,rpt):
     
    
#1: get an new rpt rpt1
#2: move it into WS super cell, rpt2
#3: copy the H(rpt2) to H(rpt1)
#3: find all dege rpts of rpt2 {rpt2}
#4: update degen of {rpt2}U{rpt1} with value n(old rpt2)+1

np.set_printoptions(precision=6,suppress=True,linewidth=150)
if __name__ == '__main__':

    # the original wannier90_hr.dat
    print ">>> build hr"

# hr, spin full
    hr = HR.from_file()

# hr_mat, no addtional spin degree
    hr_mats = hr.get_hr(0)
    nwann = hr.nwann
    norbs = hr.nwann

    dim=21;sdim=10
    degen_rpts((ndim,ndim,ndim),dtype=np.int32)

    hr_mat = np.zeros((nwann,nwann,dim,dim,dim), dtype=np.complex128) - 10000.0
    hr_mat0 = np.zeros((nwann,nwann,dim,dim,dim), dtype=np.complex128) 

# get degen in new format: easy to index
    for irpt in range(hr.nrpt):
        ix=hr.rpts[irpt,0];ixo=ixo+sdim
        iy=hr.rpts[irpt,1];iyo=iyo+sdim
        iz=hr.rpts[irpt,2];izo=izo+sdim
        degen_rpts[izo,iyo,ixo] = hr.deg_rpts[irpt]
        hr_mat[:,:,izo,iyo,ixo] = hr.mat[irpt,:,:]

# set of rpts in WS super cell
    degen_rpts_save = degen_rpts * 1

    hr_mat0 = hr_mat*1.0

# rot matrix in terms of primitive cell basis
    rot = np.zeros((3,3),np.float64)
    rotmap = np.zeros((natom),np.int32) # rot map for atoms

# -------------------------------------------------------------------------------------
# Firstly, if one of R in {Rot_R} is out of WS points set, let H({rot_R}) = 0.0
    irpt = 0
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                # skip points out of hr.rpts
                if hr_mat0[0,0,izo,iyo,ixo].real < -5000.0: 
                    continue
                rpt = np.zeros((3), dtype=np.float64)
                rpt[0]=rx*1.0;rpt[1]=ry*1.0;rpt[2]=rz*1.0
                irpt = irpt + 1
                print "--------------------irtp", irpt,  rpt,"-------------------"
                for jatom in range(natom):
                    for iatom in range(natom):
                        for isymm in range(nsymm):         
                            for ii in range(0,3):
                                for jj in range(0,3):
                                    rot[ii,jj] = np.float64(symop[isymm][0][ii][jj])
                            rot_rpt = np.dot(rot,rpt) 

                            for iptr in range(nptrans):
                                new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
                                 
                                # index of new rpt in hr_mat
                                ixt = np.int32(np.round(new_rpt[0])) 
                                iyt = np.int32(np.round(new_rpt[1])) 
                                izt = np.int32(np.round(new_rpt[2]))
                                ixn = ixt+sdim; iyn =iyt+sdim; izn = izt+sdim;
                                # determine whether new rpt out of given Huge box
                                if (abs(ixt)>10 or abs(iyt)>10 or abs(izt)>10): 
                                    print "--------------------Out of Huge Box, irtp", irpt, rpt, "-->", new_rpt, "-------------------"
                                    print "please increase dim sdim and make sure dim=2*sdim+1!"
                                    sys.exit(0)
                                # check whether the new_rpt is really new
                                if hr_mat[0,0,izn,iyn,ixn].real < -5000.0: #new rpt
                                   # find one point in orignal supercell equivalent to this point  
                                   n1,n2,n3=get_vrpt  
                                   
                                    
                               # add this point
                                    hr_mat[:,:,izn,iyn,ixn] = hr_mat0[:,:,izo,iyo,izo]
                               # update degeneracy
                                    degen_rpts[izo,iyo,ixo] = degen_rpts[izo,iyo,ixo]+1
                                    degen_rpts[izt,iyt,ixt] = degen_rpts[izo,iyo,ixo]


#----------------------------- Symmetrization ------------------------------------------------
    irpt = 0
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                if hr_mat0[0,0,izo,iyo,ixo].real < -5000.0: # points out of hr.rpts
                    continue
                rpt = np.zeros((3), dtype=np.float64)
                rpt[0]=rx*1.0;rpt[1]=ry*1.0;rpt[2]=rz*1.0
                print "--------------------irtp", irpt,  rpt,"-------------------"
                irpt = irpt + 1
                 
                for jatom in range(natom):
                    offj=2*offs[jatom] # off set of jatom in norbs 
                    norj=2*nors[jatom] # number of orbitals of jatom
                    for iatom in range(natom):
                        offi=2*offs[iatom]
                        nori=2*nors[iatom]

                        hr_ij=np.zeros((nori,norj),dtype=np.complex128)

                        for isymm in range(nsymm):         
                            for ii in range(0,3):
                                for jj in range(0,3):
                                    rot[ii,jj] = np.float64(symop[isymm][0][ii][jj])
                            rot_rpt = np.dot(rot,rpt) 

                            for iptr in range(nptrans):
                                for ii in range(0,natom):
                                    rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]

                                offjp=2*offs[rotmap[jatom]] # off set of the atom which jatom rot to in norbs
                                norjp=2*nors[rotmap[jatom]] # number of orbitals of the atom which jatom rot to 
                                offip=2*offs[rotmap[iatom]]
                                norip=2*nors[rotmap[iatom]]

                                new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
                                 
# index of new rpt in hr_mat    
                                ixt = np.int32(np.round(new_rpt[0])) 
                                iyt = np.int32(np.round(new_rpt[1])) 
                                izt = np.int32(np.round(new_rpt[2]))
 
                                ixn = ixt+sdim; iyn =iyt+sdim; izn = izt+sdim;

# determine whether new rpt     out of given rpts box
                                if abs(ixt)>10 or abs(iyt)>10 or abs(izt)>10: # for my case, ixyzt>10 may occur once in nsymm. so here should dnorm-1
                                    hr_mat_symmed[:,:,izo,iyo,ixo] = 0.0
                                    continue

                                protj = protmat[offj:offj+norj,offj:offj+norj,iptr,isymm] # collum index is of jatom, but row index is still of jatom, but with meaning of rotmap(jatom) 
                                proti = protmat[offi:offi+nori,offi:offi+nori,iptr,isymm] # protmat is blocks arranged along diagonal line. off diagonal blocks are zeros.
                                hr_ipjp= hr_mat[offip:offip+norip,offjp:offjp+norjp,izn,iyn,ixn]
                                hr_ij =  hr_ij + np.dot(np.dot(np.conj(np.transpose(proti)), hr_ipjp), protj)
# averaging
                        hr_mat_symmed[offi:offi+nori,offj:offj+norj,izo,iyo,ixo] = hr_ij/np.float64(nsymm*nptrans)

#------------------------------------------------------------------------------------------------
# only dump those with non-zero elements 
    nrpt_tmp = dim*dim*dim
    hr_mat0 = np.zeros((nrpt_tmp,norbs,norbs),dtype = np.complex128)
    rpts = np.zeros((nrpt_tmp,3),dtype = np.int32);
    dege = np.zeros((nrpt_tmp),dtype = np.int32);
    irpt = 0 # count for useful rpts
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                if abs(hr_mat_symmed[:,:,izo,iyo,ixo].real).sum() < 1.0E-9: # ommit too small 
                    continue
                dege[irpt] = 1
                rpts[irpt,0] = rx; rpts[irpt,1] = ry; rpts[irpt,2] = rz; 
                hr_mat0[irpt,:,:] = hr_mat_symmed[:,:,izo,iyo,ixo]
                irpt = irpt + 1
    nrpt = irpt
    dege_rpts = dege[0:nrpt]

#------------------------------------------------------------------------------------------------
# dump hr_mat. this should also be symmed in r space.
    print "dump new hr.dat..."
    today = datetime.date.today()
    with open('hr_symmed.dat_nsymm'+str(nsymm), 'w') as f:
        line=" Writen on "+today.ctime()+"\n"+"          "+ str(norbs) + "\n" + "        "+ str(nrpt) + "\n"
        f.write(line)
        nl = np.int32(np.ceil(nrpt/15.0))
        for l in range(nl):
            line="    "+'    '.join([str(np.int32(i)) for i in dege_rpts[l*15:(l+1)*15]])
            f.write(line)
            f.write('\n')
        for irpt in range(nrpt):
            rx = rpts[irpt,0];ry = rpts[irpt,1];rz = rpts[irpt,2]
            for jatomorb in range(norbs):
                for iatomorb in range(norbs):
                   rp =hr_mat0[irpt,iatomorb,jatomorb].real
                   ip =hr_mat0[irpt,iatomorb,jatomorb].imag
                   line="{:5d}{:5d}{:5d}{:5d}{:5d}{:18.12f}{:18.12f}\n".format(rx,ry,rz,iatomorb+1,jatomorb+1,rp,ip)	
                   f.write(line)

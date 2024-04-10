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
if __name__ == '__main__':

    # the original wannier110_hr.dat
    print ">>> copyright: Changming Yue from Department of Physics, Southern University of Science and Technology (SusTech, China). yuechangming8@gmail.com"
    print ">>> Loading wannier90_hr.dat ..." 

    # max hopping distance in terms of lattice vector -sdim to sdim
    # allowed by Rot( all Rpts ). Please reset it if out of huge box
    sdim = 16

    # hdims of rpts
    ndim = 2*sdim+1

    hr = HR.from_file_large(ndim, sdim,'wannier90_hr.dat')

# hr_mat(norbs,norbs,ndim,ndim) j,i,rz,ry,rx  <i|H(rx,ry,rz)|j>
    hr_mat = hr.get_hr_large()  # forced zeros version 

# backuped for checking whether {rot_R} is out of wannier90_hr.dat's rpts
    hr_mat0= hr_mat * 1.0 # times 1.0 in case of quotation
    hr_mats= hr_mat * 1.0 # times 1.0 in case of quotation

    print "ispinor", ispinor

#-------------- #-------------- #-------------- #-------------- #-------------- #--------------
#-------------- #-------------- #-------------- #-------------- #-------------- #--------------

    print "Make Hemitrization and Time reversal symmetry ...." 
    if ispinor:  # for spinor case, norbs=nband*2. 
        print "------------------spinfull--------------------------------"
        print "####### It is assumed up up dn dn wannier90_hr.dat ########"
        print "Temporarily up up dn dn will be changed up dn up dn."
        print "----------------------------------------------------------"
        # pauli matrix (Sigma_y)^T * 1.0j, Sigma_y * 1.0j
        syl = np.zeros((2,2),dtype = np.complex128)
        syr = np.zeros((2,2),dtype = np.complex128)
        syl[0,1]=-1.0;syl[1,0]=1.0;syr[0,1]=1.0;syr[1,0]=-1.0;
        umata = np.zeros((norbs,norbs), dtype = np.complex128)
        umatl = np.zeros((norbs,norbs), dtype = np.complex128)
        umatb = np.zeros((norbs,norbs), dtype = np.complex128)
        umatr = np.zeros((norbs,norbs), dtype = np.complex128)
        for iorbs in range(nband):
            umata[2*iorbs:2*iorbs+2,2*iorbs:2*iorbs+2] = syl
            umatb[2*iorbs:2*iorbs+2,2*iorbs:2*iorbs+2] = syr
        # up up
        umatl[0:nband,0:nband] = umata[0:norbs:2,0:norbs:2]
        # up dn
        umatl[0:nband,nband:norbs] = umata[0:norbs:2,1:norbs:2]
        # dn up
        umatl[nband:norbs,0:nband] = umata[1:norbs:2,0:norbs:2]
        # dn dn
        umatl[nband:norbs,nband:norbs] = umata[1:norbs:2,1:norbs:2]
    
        # up up
        umatr[0:nband,0:nband] = umatb[0:norbs:2,0:norbs:2]
        # up dn
        umatr[0:nband,nband:norbs] = umatb[0:norbs:2,1:norbs:2]
        # dn up
        umatr[nband:norbs,0:nband] = umatb[1:norbs:2,0:norbs:2]
        # dn dn
        umatr[nband:norbs,nband:norbs] = umatb[1:norbs:2,1:norbs:2]
    
        # -------------------------------- Hermitrization -------------------------------------
        print "Hermitization restore ..."
        irpt = 0
        for rx in range(-sdim,sdim+1):
            #print "rx ", rx
            ixo = rx+sdim
            for ry in range(-sdim,sdim+1):
                iyo = ry+sdim
                for rz in range(-sdim,sdim+1):
                    izo = rz+sdim

                    if hr_mats[0,0,izo,iyo,ixo].real < -5000.0: 
                        continue
    
                    miz = -rz + sdim
                    miy = -ry + sdim
                    mix = -rx + sdim
    
                    if hr_mats[0,0,miz,miy,mix].real < -5000.0: 
                        print "Error! rpts not symmetric in original r-mesh!!"
                        sys.exit(0)
                    hr_mat[:,:,izo,iyo,ixo] = (hr_mat0[:,:,izo,iyo,ixo] + np.conj(np.transpose(hr_mat0[:,:,miz,miy,mix])))/2.0
    
                    # qswu think that for rpts with dege>1, actual hopping should be smaller. 
                    hr_mat[:,:,izo,iyo,ixo] = hr_mat[:,:,izo,iyo,ixo] / hr.deg_rpt[irpt]
                    hr.deg_rpt[irpt]=1
                    irpt = irpt + 1
    
        # ---------------------------- Time reversal -------------------------------------------
        print "time reversal restore ..."
        for rx in range(-sdim,sdim+1):
            #print "rx ", rx
            ixo = rx+sdim
            for ry in range(-sdim,sdim+1):
                iyo = ry+sdim
                for rz in range(-sdim,sdim+1):
                    izo = rz+sdim
                    if hr_mats[0,0,izo,iyo,ixo].real < -5000.0: 
                        continue
                    old = hr_mat[:,:,izo,iyo,ixo] * 1.0
                    new = np.conj(hr_mat[:,:,izo,iyo,ixo]) * 1.0
                    new = np.dot(np.dot(umatl,new),umatr) 
                    hr_mat[:,:,izo,iyo,ixo] = (old + new)/2.0
    
                    # adjust data structure
                    hr_mat0[0:norbs:2,0:norbs:2,izo,iyo,ixo] = hr_mat[0:nband,0:nband,izo,iyo,ixo]
                    # up dn
                    hr_mat0[0:norbs:2,1:norbs:2,izo,iyo,ixo] = hr_mat[0:nband,nband:norbs,izo,iyo,ixo]
                    # dn up
                    hr_mat0[1:norbs:2,0:norbs:2,izo,iyo,ixo] = hr_mat[nband:norbs,0:nband,izo,iyo,ixo]
                    # dn dn
                    hr_mat0[1:norbs:2,1:norbs:2,izo,iyo,ixo] = hr_mat[nband:norbs,nband:norbs,izo,iyo,ixo]
        hr_mat = hr_mat0 * 1.0

    else:
        print "------------------spinless--------------------------------"
	irpt = 0
        for rx in range(-sdim,sdim+1):
            ixo = rx+sdim
            for ry in range(-sdim,sdim+1):
                iyo = ry+sdim
                for rz in range(-sdim,sdim+1):
                    izo = rz+sdim
                    miz = -rz + sdim
                    miy = -ry + sdim
                    mix = -rx + sdim
                    hr_mat[:,:,izo,iyo,ixo] = (hr_mat0[:,:,izo,iyo,ixo] + np.conj(np.transpose(hr_mat0[:,:,miz,miy,mix])))/2.0
                    # imag part should be forced to zeros for spin-less real-space H on real-orbitals basis
                    hr_mat[:,:,izo,iyo,ixo] = hr_mat[:,:,izo,iyo,ixo].real
                    # qswu think that for rpts with dege>1, actual hopping should be smaller. 
                    if hr_mats[0,0,izo,iyo,ixo].real > -5000.0: 
                        hr_mat[:,:,izo,iyo,ixo] = hr_mat[:,:,izo,iyo,ixo] / hr.deg_rpt[irpt]
			hr.deg_rpt[irpt] = 1
			irpt = irpt + 1

#-------------- #-------------- #-------------- #-------------- #-------------- #--------------
#-------------- #-------------- #-------------- #-------------- #-------------- #--------------

    hr_mat0 = hr_mat * 1.0 # times 1.0 in case of quotation
    hr_mats = []

# rot matrix in terms of primitive cell basis
    rot = np.zeros((3,3),np.float64)
    rotmap = np.zeros((natom_wan),np.int32) # rot map for atoms
    rdege = np.zeros((ndim,ndim,ndim),dtype=np.int32) + 1

    print ">>> Performing symmerizations, for details please refer symm.log"
    tenpct = np.int32(np.linspace(0,hr.nrpt,12))

# flags for new added block or not: =1 for yes
    flags = np.zeros((natom_wan,natom_wan,ndim,ndim,ndim),dtype=np.int32) 

# -------------------------------------------------------------------------------------
    print ">>> Step 1, find new rpts block by rotate wannier90_hr.dat's rpt blocks"
    irpt = 0
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                # rpts out of wannier90_hr.dat's rpts. we get new rpt block by existed ones.
                if hr_mat0[0,0,izo,iyo,ixo].real < -5000.0: 
                    continue
                rpt = np.zeros((3), dtype=np.float64)
                rpt[0]=rx*1.0;rpt[1]=ry*1.0;rpt[2]=rz*1.0
                if (irpt in tenpct): print "{:12.6f}".format(1.0*irpt/hr.nrpt*100),"%..."
                irpt = irpt + 1
                 
                for jatom in range(natom_wan):
                    for iatom in range(natom_wan):
                        for isymm in range(nsymm):         
                            rot = np.float64(symop[isymm][0])
                            rot_rpt = np.dot(rot,rpt) 

                            for iptr in range(nptrans):
                                for ii in range(0,natom_wan):
                                    rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]

                                new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
                                 
                                ixt = np.int32(np.round(new_rpt[0])) 
                                iyt = np.int32(np.round(new_rpt[1])) 
                                izt = np.int32(np.round(new_rpt[2]))
 
                                ixn = ixt+sdim; iyn =iyt+sdim; izn = izt+sdim;

                                if abs(ixt)>sdim or abs(iyt)>sdim or abs(izt)>sdim: 
				    print "--------------------Out of Huge Box, irtp", irpt, rpt, "-->", new_rpt, "-------------------"
				    print "Error! please open symmhrspinless.py and make a bigger sdim!!!!!! "
				    sys.exit(0)
                                else: 
                                    if hr_mat0[0,0,izn,iyn,ixn].real < -5000.0:
				        # we find new rpt blocks 
				        flags[rotmap[iatom],rotmap[jatom],izn,iyn,ixn] = 1

    print "I find the num of new blocks is:", sum(sum(sum(sum(sum(flags))))), "These new rpt blocks are\
    listed as below. Meanwhile, hoppings of these new blocks will obtained from wannier90_hr.dat."

# -------------------------------------------------------------------------------------
    # prefix of dims, for spinor case spin-orbital is doubled
    ndm = 1
    if ispinor: ndm = 2

# get hoppings of new rpt block from wannier90_hr.dat existed (stored in hr_mat0)
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                rpt = np.zeros((3), dtype=np.float64)
                rpt[0]=rx*1.0;rpt[1]=ry*1.0;rpt[2]=rz*1.0
                for jatom in range(natom_wan):
                    offj=ndm*offs[jatom] # off set of jatom in norbs 
                    norj=ndm*nors[jatom] # number of orbitals of jatom
                    for iatom in range(natom_wan):
                        offi=ndm*offs[iatom]
                        nori=ndm*nors[iatom]
                        if flags[iatom,jatom,izo,iyo,ixo]==1:
     		           print "--new rpt block--", rx,ry,rz,iatom,jatom 
			   cnt = 0 
			   lexit = 0 # once found and got, exit
                           for isymm in range(nsymm):         
			      if lexit==1: break
                              rot = np.float64(symop[isymm][0])
                              rot_rpt = np.dot(rot,rpt) 

                              for iptr in range(nptrans):
			         if lexit==1: break
                                 for ii in range(0,natom_wan):
                                     rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]

                                 offjp=ndm*offs[rotmap[jatom]] # off set of the atom which jatom rot to in norbs
                                 norjp=ndm*nors[rotmap[jatom]] # number of orbitals of the atom which jatom rot to 
                                 offip=ndm*offs[rotmap[iatom]]
                                 norip=ndm*nors[rotmap[iatom]]

                                 new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
                                 ixt = np.int32(np.round(new_rpt[0])) 
                                 iyt = np.int32(np.round(new_rpt[1])) 
                                 izt = np.int32(np.round(new_rpt[2]))
                                 ixn = ixt+sdim; iyn =iyt+sdim; izn = izt+sdim;

                                 if abs(ixt)>sdim or abs(iyt)>sdim or abs(izt)>sdim: 
			             print "Error: out of huge box!!!"
			             sys.exit(0)

				 # maybe it is a new rpt. There must be a block in wannier90_hr.dat rotated 
				 # to this block. if cnt==nsymm*nptran and lexit==0 error!!
                                 if hr_mat0[0,0,izn,iyn,ixn].real < -5000.0: 
			             cnt = cnt+1
                                     continue

                                 protj = protmat[offj:offj+norj,offj:offj+norj,iptr,isymm] # 
                                 proti = protmat[offi:offi+nori,offi:offi+nori,iptr,isymm] # 
                                 hr_ipjp= hr_mat[offip:offip+norip,offjp:offjp+norjp,izn,iyn,ixn]
                                 hr_ij =  np.dot(np.dot(np.conj(np.transpose(proti)), hr_ipjp), protj)
                                 hr_mat[offi:offi+nori,offj:offj+norj,izo,iyo,ixo] = hr_ij
			 	 cnt = cnt + 1
			  	 lexit = 1
                           if lexit==0: 
			         print "Error: rot not close!!"
			         print "cnt", cnt
				 sys.exit(0)

# ----------------------------------------------------------------------------------------------------
# loop over new generated rpts. set new rpts' block to be zero if that block is not swept  
# since not all block in new rpt has correpondence to wannier90_hr.dat's. 
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
		# a new rpts
                s = flags[:,:,izo,iyo,ixo]
		sm = sum(sum(s))
                if sm > 0: # this rpt is a new rpt
		     for i in range(norbs):
		         for j in range(norbs): 
			     a  = hr_mat[i,j,izo,iyo,ixo].real
			     if a < -5000.0: hr_mat[i,j,izo,iyo,ixo] = 0.0	

# store hr_mat again to hr_mat0
    hr_mat0 = hr_mat * 1.0

#----------------------------- Symmetrization ------------------------------------------------
    print "Step 2, averaing to symmetrize"
    hr_mat_symmed = hr_mat * 0.0
    irpt = 0
    for rx in range(-sdim,sdim+1):
        ixo = rx+sdim
        for ry in range(-sdim,sdim+1):
            iyo = ry+sdim
            for rz in range(-sdim,sdim+1):
                izo = rz+sdim
                # points out wannier90_hr.dat's rpts + new rpts
                if hr_mat0[0,0,izo,iyo,ixo].real < -5000.0: 
                    continue
                rpt = np.zeros((3), dtype=np.float64)
                rpt[0]=rx*1.0;rpt[1]=ry*1.0;rpt[2]=rz*1.0
                #print>>f, "--------------------irtp", irpt,  rpt,"-------------------"
                if (irpt in tenpct): print "{:12.6f}".format(1.0*irpt/hr.nrpt*100),"%..."
                irpt = irpt + 1
                 
                for jatom in range(natom_wan):
                    offj=ndm*offs[jatom] # off set of jatom in norbs 
                    norj=ndm*nors[jatom] # number of orbitals of jatom
                    for iatom in range(natom_wan):
                        offi=ndm*offs[iatom]
                        nori=ndm*nors[iatom]

                        hr_ij=np.zeros((nori,norj),dtype=np.complex128)
			cnt = 0
			chug = 0
                        for isymm in range(nsymm):         
                            rot = np.float64(symop[isymm][0])
                            rot_rpt = np.dot(rot,rpt) 

                            for iptr in range(nptrans):
                                for ii in range(0,natom_wan):
                                    rotmap[ii]=wann_atom_rotmap[isymm][iptr][ii]

                                offjp=ndm*offs[rotmap[jatom]] # off set of the atom which jatom rot to in norbs
                                norjp=ndm*nors[rotmap[jatom]] # number of orbitals of the atom which jatom rot to 
                                offip=ndm*offs[rotmap[iatom]]
                                norip=ndm*nors[rotmap[iatom]]

                                new_rpt = rot_rpt + vec_shift[isymm,iptr,jatom] - vec_shift[isymm,iptr,iatom] 
                                ixt = np.int32(np.round(new_rpt[0])) 
                                iyt = np.int32(np.round(new_rpt[1])) 
                                izt = np.int32(np.round(new_rpt[2]))
                                ixn = ixt+sdim; iyn =iyt+sdim; izn = izt+sdim;
                 
		                # blocks in new rpts may rot out of huge box or out of wann's rpts+new rpts 
				# set this block zero ( alread set in init of hr_symmed.dat )
                                if abs(ixt)>sdim or abs(iyt)>sdim or abs(izt)>sdim: 
				    cnt = cnt + 1
				    chug = chug + 1
                                    continue
                                if hr_mat0[0,0,izn,iyn,ixn].real < -5000.0: 
				    cnt = cnt + 1
                                    continue

                                protj = protmat[offj:offj+norj,offj:offj+norj,iptr,isymm] 
                                proti = protmat[offi:offi+nori,offi:offi+nori,iptr,isymm] 
                                hr_ipjp= hr_mat[offip:offip+norip,offjp:offjp+norjp,izn,iyn,ixn]
                                hr_ij =  hr_ij + np.dot(np.dot(np.conj(np.transpose(proti)), hr_ipjp), protj)

                        tmp = hr_mat[offi:offi+nori,offj:offj+norj,izo,iyo,ixo]
			if cnt > 1 and cnt < nsymm*nptrans and sum(sum(abs(tmp)))>0.0: 
			   print "Error!!!", cnt, chug, sum(sum(abs(tmp)))
			   sys.exit(0)
                           
                        hr_mat_symmed[offi:offi+nori,offj:offj+norj,izo,iyo,ixo] = hr_ij/np.float64(nsymm*nptrans)

#------------------------------------------------------------------------------------------------
# only dump those with non-zero elements 
    nrpt_tmp = ndim*ndim*ndim
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
                dege[irpt] = rdege[izo,iyo,ixo]
                rpts[irpt,0] = rx; rpts[irpt,1] = ry; rpts[irpt,2] = rz; 
                hr_mat0[irpt,:,:] = hr_mat_symmed[:,:,izo,iyo,ixo]
                irpt = irpt + 1
    nrpt = irpt
    dege_rpts = dege[0:nrpt]

    # deallocate
    hr_mat_symmed = []
    
    # for spinor case change up dn up dn to up up dn dn 
    if ispinor: 
        hr_mat = hr_mat0 * 0.0
        for irpt in range(nrpt):
            hr_mat[irpt,0:nband,0:nband]     = hr_mat0[irpt,0:norbs:2,0:norbs:2]
            hr_mat[irpt,0:nband,nband:norbs] = hr_mat0[irpt,0:norbs:2,1:norbs:2]
            hr_mat[irpt,nband:norbs,0:nband] = hr_mat0[irpt,1:norbs:2,0:norbs:2]
            hr_mat[irpt,nband:norbs,nband:norbs] = hr_mat0[irpt,1:norbs:2,1:norbs:2]
    else:
        hr_mat = hr_mat0 * 1.0

#------------------------------------------------------------------------------------------------
# dump hr_mat. this should also be symmed in r space.
    print "Step 3, dump new hr.dat..."
    today = datetime.date.today()
    dt = str(datetime.datetime.now())
    with open('wannier90_hr.dat_nsymm'+str(nsymm), 'w') as f:
        line="Symmed by ChangmingYue,yuechangming8@gamil.com, on "+dt+"\n"+"          "+ str(norbs) + "\n" + "        "+ str(nrpt) + "\n"
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
                   rp =hr_mat[irpt,iatomorb,jatomorb].real
                   ip =hr_mat[irpt,iatomorb,jatomorb].imag
                   line="{:5d}{:5d}{:5d}{:5d}{:5d}{:20.14f}{:20.14f}\n".format(rx,ry,rz,iatomorb+1,jatomorb+1,rp,ip)	
                   f.write(line)

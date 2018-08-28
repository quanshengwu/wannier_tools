#!/usr/bin/env python
'''
Developer: Changming Yue, yuechangming8@gmail.com
perform symmetrization for tight-binding Hamiltonian abtained by 
wannier90 *in momentum space*. For equivalent atoms, local axices
should be set in wannier90.win file. 
'''
import datetime,sys,io
import numpy as np

from wanntb.kvec import SymKVec, UniKVec
from wanntb.tran import fourier_hr2h1k, fourier_hr2hk, tran_op
from wanntb.io import write_band
from progressbar import *

from numpy.linalg import eigh 
from wanntb.ham import HR
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# define types and size
Amat=np.zeros((3,3),dtype=np.float64)
# ------------------------------Change it according your case------------------------
# primitive cell basis 

Amat[:,0]= [   -3.8889403198684809,   2.5564011223049747,   4.1324340392468253] 
Amat[:,1]= [    3.8889403198600525,  -2.5564011223045635,   4.1324340392455445]
Amat[:,2]= [    3.8889403198630697,   2.5564011222967848,  -4.1324340392466352]
 


# ------------------------------Change it according your case------------------------
# high symmetry K-point in 1st BZ.
hsymkpt=[[0.00,	 0.00,	 0.000],# W
         [0.75,	-0.25,	-0.250],# L
         [0.50,	 0.00,	 0.000],# G
         [0.00,	 0.00,	 0.000],# X
         [0.50,	 0.50,	-0.500],# X
         [0.50,	 0.00,	-0.500],# X
         [0.00,	 0.00,	 0.000]]# W

# number of kpoints on a k-path, say from G-X
nkpt_path = 50
# --------------------------End-Change it according your case------------------------
 
# more precise than kbase by direct calc from lattice vectors
kvec=np.zeros((3,3),dtype=np.float64)
vol = np.dot(Amat[:,0],np.cross(Amat[:,1],Amat[:,2]))
kvec[:,0] = np.cross(Amat[:,1],Amat[:,2])/vol
kvec[:,1] = np.cross(Amat[:,2],Amat[:,0])/vol
kvec[:,2] = np.cross(Amat[:,0],Amat[:,1])/vol

print "please use kvec instead of kbase in generating kpts to improve precision!"
print "kvec.T:\n", kvec.T*np.pi*2

np.set_printoptions(precision=6,suppress=True,linewidth=150)
if __name__ == '__main__':

    # the original wannier90_hr.dat
    if rank==0:print ">>> build hr"

    # hr, no spin freedom
    hr = HR.from_file('hr_symmed.dat_nsymm8')
    #hr = HR.from_file('wannier90_hr.dat')

#  dimensions
    nband = hr.nwann     
    norbs = nband
    print "my rank:", rank

    # hr_mat with spin freedom, all up terms makes block of wannier90_hr.dat
    # spin-less
    hr_mat0 = hr.get_hr(0)

# define a symmetry k-path
    symk = SymKVec(np.float64(kvec), np.float64(hsymkpt))
    symk.from_hsymkpt(nkpt_path)
    symk.get_klen()

    hk = np.zeros((symk.nkpt,norbs,norbs), dtype=np.complex128)

# loop over all kpoints, FT with gauge factor
    if rank==0:print "FT H(R)->H(k) with gauge factor exp(i k*tau) "
    if rank==0:pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=symk.nkpt).start()
    for ikpt in range(rank,symk.nkpt,size):
        rkpt = symk.kvec[ikpt]
        hk0  = fourier_hr2h1k(norbs, rkpt, hr.nrpt, hr.rpts, hr.deg_rpt, hr_mat0)

        hk[ikpt,:,:] = hk0
        if rank==0:pbar.update(ikpt)
    if rank==0:pbar.finish()
    hk = comm.allreduce(hk)
    comm.barrier()

# dump energy band without soc
    if rank==0:print ">>> get band structure without soc"
    eigval = np.zeros((symk.nkpt, norbs),dtype=np.float64)
    eigvec = np.zeros((symk.nkpt, norbs, norbs),dtype=np.complex128)
    if rank==0:pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=symk.nkpt).start()
    for i in range(rank,symk.nkpt,size):
        eigval[i,:], eigvec[i,:,:] = eigh(hk[i,:,:])
        if rank==0:pbar.update(i)
    if rank==0:pbar.finish()
    eigval = comm.allreduce(eigval)
    comm.barrier()


# dump wann soc  band
    if rank==0:write_band(symk.klen, eigval, 'wann.band.dat')

#!/usr/bin/env python

import numpy as np
from numpy.linalg import eigh
from tran import fourier_hr2hk
import sys

class HR():
    """ 
    The Wannier tight-binding (TB) Hamiltonian in real space.

    """

    def __init__(self, nwann, nrpt, irpt0, rpts, deg_rpt, hr):
        self.nwann=nwann
        self.nrpt = nrpt
        self.irpt0 = irpt0 
        self.deg_rpt = deg_rpt
        self.rpts = rpts
        self.hr = hr

    @staticmethod
    def from_file(fname='wannier90_hr.dat'):

        """
        Generate a TB Hamiltonian from Wannier90 output file "case_hr.dat".

        Args:
            fname (string): the file that contains the Wannier90 output file 
            "case_hr.dat"
        """

        try: 
            with open(fname, 'r') as f:
                # skip the header (1 line)
                f.readline()
                nwann = int(f.readline().strip()) 
                nrpt = int(f.readline().strip())
                # the degeneracy of R points
                nline = np.int32(np.ceil(nrpt/15.0))
                print "nline", nline
                tmp = []
                for i in range(nline):
                    tmp.extend(f.readline().strip().split())   
                tmp = [np.int(item) for item in tmp]
                deg_rpt = np.array(tmp, dtype=np.int)
                # read hr for each r-point
                rpts = np.zeros((nrpt,3), dtype=np.int) 
                hr = np.zeros((nrpt, nwann, nwann), dtype=np.complex128)
                for i in range(nrpt):
                    print "i", i
                    for j in range(nwann):
                        for k in range(nwann):
                            rx, ry, rz, hr_i, hr_j, hr_real, hr_imag = f.readline().strip().split() 
                            rpts[i,:] = int(rx), int(ry), int(rz) 
                            if int(rx) == 0 and int(ry) == 0 and int(rz) == 0:
                                irpt0 = i  
                            hr[i,k,j] = np.float64(hr_real) + np.float64(hr_imag) * 1j
                # construct the HR instance
                return HR(nwann, nrpt, irpt0, rpts, deg_rpt, hr)

        except IOError:
            print "File:" + "\"" + fname + "\"" +  " doesn't exist!"

    @staticmethod
    def from_file_large(dim,sdim,fname='wannier90_hr.dat'):

        """
        Generate a TB Hamiltonian from Wannier90 output file "case_hr.dat".
        We use index directly as rpt up to a const vector
        Args:
            fname (string): the file that contains the Wannier90 output file 
            "case_hr.dat"
        """
        print "in from_file_large dim, sdim", dim, sdim
        try: 
            #fh_info = open('deg_rpt','w')
            with open(fname, 'r') as f:
                # skip the header (1 line)
                f.readline()
                nwann = int(f.readline().strip()) 
                nrpt = int(f.readline().strip())
                print "nrpt nwann", nrpt, nwann
                if nrpt > dim*dim*dim:
                      print "Too much rpts in your hr file!"
                      sys.exit(0)
                # the degeneracy of R points
                nline = np.int32(np.ceil(nrpt/15.0))
                tmp = []
                for i in range(nline):
                    tmp.extend(f.readline().strip().split())   
                tmp = [np.int(item) for item in tmp]
                deg_rpt = np.array(tmp, dtype=np.int)
                rpts = np.zeros((nrpt,3), dtype=np.int) 
                hr = np.zeros((nwann,nwann,dim,dim,dim), dtype=np.complex128)
                hr = hr - 10000.0
                tenpct = np.int32(np.linspace(0,nrpt,12))
                print "loading hamr..."
                for i in range(nrpt):
                    if (i in tenpct): print "{:12.6f}".format(1.0*i/nrpt*100),"%..."
                    for j in range(nwann):
                        for k in range(nwann):
                            rx, ry, rz, hr_i, hr_j, hr_real, hr_imag = f.readline().strip().split() 
                            rpts[i,:] = int(rx), int(ry), int(rz) 
                            if abs(int(rx)) > sdim or abs(int(rx)) > sdim or  abs(int(rx)) > sdim :
                                print "rpts index out of range of [-sdim:sdim]", rx, ry, rz, "sdim", sdim
                                sys.exit(0)

                            idx = np.int(rx)+sdim; idy = np.int(ry)+sdim; idz=np.int(rz)+sdim
# remember when dumpy print idx idy, idz k, j
                            hr[k,j,idz,idy,idx]=np.float64(hr_real) + np.float64(hr_imag) * 1j
                    #print >>fh_info, "rx ry rx", rx, ry, rz, "dege", deg_rpt[i]
                # construct the HR instance
                print "loading done"
                return HR(nwann, nrpt, 0, rpts, deg_rpt, hr)

        except IOError:
            print "File:" + "\"" + fname + "\"" +  " doesn't exist!"

    @staticmethod
    def copy_hr(other):

        """
        copy instance of HR
        """

        return HR( other.nwann, 
                   other.nrpt, 
                   other.irpt0, 
                   np.copy(other.rpts), 
                   np.copy(other.deg_rpt), 
                   np.copy(other.hr) )

    def get_hr0(self, ispin=False):
    
        """
        return the on-site term 
        """

        if ispin:
            norbs = 2 * self.nwann
            hr0_spin = np.zeros((norbs, norbs), dtype=np.complex128) 
            hr0_spin[0:norbs:2, 0:norbs:2] = self.hr[self.irpt0,:,:]
            hr0_spin[1:norbs:2, 1:norbs:2] = self.hr[self.irpt0,:,:]
            return hr0_spin
        else:
            return self.hr[self.irpt0,:,:]
            
    def get_hr(self, ispin):
        
        """
        return hr
        """

        # with spin, spin order: up dn up dn up dn
        if ispin==1:
            #print "for spin-full and time reversal invariant case, dn block is conj of up block"
            norbs = 2 * self.nwann 
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:norbs:2, 0:norbs:2] = self.hr
            hr_spin[:, 1:norbs:2, 1:norbs:2] = self.hr
            return hr_spin
        # with spin, spin order: up up up dn dn dn
        elif ispin==2:
            norbs = 2 * self.nwann 
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:self.nwann, 0:self.nwann] = self.hr
            hr_spin[:, self.nwann:norbs, self.nwann:norbs] = self.hr
            return hr_spin
        # without spin
        else:
            return self.hr[:,:,:]

    def get_hr_large(self):
        
        """
        return hr
        """

        return self.hr[:,:,:,:,:]

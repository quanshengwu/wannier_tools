#!/usr/bin/env python

import numpy as np

class KVec():

    """
    define K points in BZ, high symmetry line or uniform grid
    """

    def __init__(self, kpt_type='uni', kbase=None, nkpt=None, kvec=None):
        self.nkpt = nkpt
        self.kbase = np.array(kbase, dtype=np.float64)
        self.kvec = np.array(kvec, dtype=np.float64)
        self.kpt_type = kpt_type

    def set_base(self,kbase):
        self.kbase = np.array(kbase, dtype=np.float64)

    def kvec_from_file(self,fname):
        tmp = []
        try: 
            with open(fname, 'r') as f:
                for line in f:
                    line = line.strip().split()
                    if line != []:
                        tmp.append(line)
                self.kvec = np.array(tmp, dtype=np.float64)
                self.nkpt = len(tmp)
        except IOError:
            print "File" + "\"" + fname + "\"" + "doesn't exists!"

 
class SymKVec(KVec):
    
    """
    K points in high symmetry line 
    """

    def __init__(self, kbase=None, hsymkpt=None, klen=None):
        self.klen = np.array(klen, dtype=np.float64)
        self.hsymkpt = np.array(hsymkpt, dtype=np.float64)
        KVec.__init__(self, 'sym', kbase)

    def get_klen(self):
        self.klen = np.zeros(self.nkpt, dtype=np.float64)          
        self.klen[0] = 0.0
        prev_kpt = self.kvec[0]
        
        for i in range(1,self.nkpt):
            curr_kpt = self.kvec[i,:]
            tmp_kpt = curr_kpt - prev_kpt
            kx = np.dot(tmp_kpt, self.kbase[:,0])
            ky = np.dot(tmp_kpt, self.kbase[:,1])
            kz = np.dot(tmp_kpt, self.kbase[:,2])
            self.klen[i] = self.klen[i-1] + np.sqrt(np.dot((kx,ky,kz), (kx,ky,kz)))
            prev_kpt = curr_kpt

    def from_hsymkpt(self, nkpt_per_path=20):
        self.nkpt = nkpt_per_path * (len(self.hsymkpt)-1)
        self.kvec = np.zeros((self.nkpt,3), dtype=np.float64)
        for i in range(1,len(self.hsymkpt)):
            kpt_prev = self.hsymkpt[i-1,:]
            kpt_curr = self.hsymkpt[i  ,:]
            for j in range(nkpt_per_path):
                ikpt = (i-1)*nkpt_per_path + j 
                self.kvec[ikpt,:] = np.float64(j)/np.float64(nkpt_per_path-1) * (kpt_curr - kpt_prev) + kpt_prev     


    def from_myhsymkpt(self): # only hsymkpt for soc fit
        self.nkpt = len(self.hsymkpt)
        self.kvec = np.zeros((self.nkpt,3), dtype=np.float64)
        for i in range(len(self.hsymkpt)):
            self.kvec[i,:] = self.hsymkpt[i,:]

    def from_hsymkpt_uni(self, step):
        kvec = []
        self.hsym_dis = np.zeros(len(self.hsymkpt),dtype=np.float64)
        self.hsym_dis[0]=0.0
        for i in range(0,len(self.hsymkpt)-1):
            kpt_prev = self.hsymkpt[i,:]
            kpt_curr = self.hsymkpt[i+1,:]
            tmp = np.dot(self.kbase.transpose(), kpt_curr-kpt_prev)
            dis = np.sqrt(np.dot(tmp,tmp))
            self.hsym_dis[i+1]=self.hsym_dis[i] + dis
            pts = np.arange(0,dis,step) / dis
            for ipt in pts:
                kvec.append(ipt * (kpt_curr-kpt_prev) + kpt_prev)
        self.kvec = np.array(kvec,dtype=np.float64)
        self.nkpt = len(self.kvec)
    

class UniKVec(KVec):

    """
    K points in uniform grid
    """

    def __init__(self, grid=None):
        self.grid = grid
        KVec.__init__(self, 'uni')

    def from_grid(self):
        delta=0.001
        nx,ny,nz = self.grid
        self.nkpt = nx * ny * nz
        self.kvec = np.zeros((self.nkpt, 3), dtype=np.float64)
        ikpt = 0
        for i in range(nx):
            if nx==1:
                kx=0.0
            else:
                kx = float(i)/float(nx)  
            for j in range(ny):
                if ny==1:
                    ky=0.0
                else:
                    ky = float(j)/float(ny)  
                for k in range(nz):
                    if nz==1:
                        kz=0.0
                    else:
                        kz = float(k)/float(nz)  
                    ikpt = ikpt + 1
                    self.kvec[ikpt-1,:] = kx+delta,ky+delta,kz+delta

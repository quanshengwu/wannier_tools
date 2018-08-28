#!/usr/bin/env python
'''
Developer: Changming Yue, yuechangming10@gmail.com
perform symmetrization for tight-binding Hamiltonian abtained by 
wannier110. For equivalent atoms, local axices should be set in
wannier110.win file. 
'''
import datetime,sys,io
import numpy as np
from wanntb.ham import HR

np.set_printoptions(precision=6,suppress=True,linewidth=150)
if __name__ == '__main__':
    hr = HR.from_file(fname='hr_symmed.dat_nsymm4')
    nband = hr.nwann/2
    nwann = hr.nwann

# for check {rot_R} and dump
    hr_mat0= hr.get_hr(0)
    hr_mat = hr_mat0 * 0.0 # times 1.0 in case of quotation of hr_mat

# ---------------------------- Time reversal -------------------------------------------
    for irpt in range(hr.nrpt):
        hr_mat[irpt,0:nband,0:nband]     = hr_mat0[irpt,0:nwann:2,0:nwann:2]
        hr_mat[irpt,0:nband,nband:nwann] = hr_mat0[irpt,0:nwann:2,1:nwann:2]
        hr_mat[irpt,nband:nwann,0:nband] = hr_mat0[irpt,1:nwann:2,0:nwann:2]
        hr_mat[irpt,nband:nwann,nband:nwann] = hr_mat0[irpt,1:nwann:2,1:nwann:2]

#------------------------------------------------------------------------------------------------
# dump hr_mat
    print "dump hr uudd ..."
    ctime = datetime.datetime.now()
    with open('wann_symm_uudd_hr.dat', 'w') as f:
        line=" Writen on "+str(ctime)+"\n"+"          "+ str(nwann) + "\n" + "        "+ str(hr.nrpt) + "\n"
        f.write(line)
        nl = np.int32(np.ceil(hr.nrpt/15.0))
        for l in range(nl):
            line="    "+'    '.join([str(np.int32(i)) for i in hr.deg_rpt[l*15:(l+1)*15]])
            f.write(line)
            f.write('\n')

# for WS points with all Hr-elements zero, we do not dump
        for irpt in range(hr.nrpt):
             rx = hr.rpts[irpt,0];ry = hr.rpts[irpt,1];rz = hr.rpts[irpt,2]
             for jatomorb in range(nwann):
                 for iatomorb in range(nwann):
                    rp =hr_mat[irpt,iatomorb,jatomorb].real
                    ip =hr_mat[irpt,iatomorb,jatomorb].imag
                    line="{:5d}{:5d}{:5d}{:5d}{:5d}{:20.12f}{:20.12f}\n".format(rx,ry,rz,iatomorb+1,jatomorb+1,rp,ip)	
                    f.write(line)

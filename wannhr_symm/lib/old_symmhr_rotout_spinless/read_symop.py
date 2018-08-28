#!/usr/bin/env python
import numpy as np
from input_paras import *
symop=[]

######################################################################
# symmop.dat is copyied from vasp's OUTCAR, and it is formmated as
# 
# irot  :   1
# --------------------------------------------------------------------
# isymop:   1   0   0
#           0   1   0
#           0   0   1

# gtrans:     0.0000000     0.0000000     0.0000000

# ptrans:     0.0000000     0.0000000     0.0000000
# rotmap:
# (   1->   1)  (   2->   2)  (   3->   3)  (   4->   4)  (   5->   5) 
# (   6->   6)  (   7->   7) 
#
#
# irot  :   2
# .
# .
######################################################################

f=open('symmop.dat','r')
for i in range(nsymm):
    f.readline() 
    f.readline() 
    rot=[]
    line=f.readline()
    a=line.split()[1:4]
    rot.append(a)
    line=f.readline()
    rot.append(line.split())
    line=f.readline()
    rot.append(line.split())

    f.readline() 
    
    line=f.readline()
    gtrans = line.split()[1:4]  
    f.readline()

    rotmap=[]
    ptrans=[]
    for nptr in range(nptrans):
        line=f.readline()
        ptrans.append(line.split()[1:4])
        f.readline()

        mp = []
# number of lines for map of atoms
        nl = np.int32(np.ceil((natom_all-0.01)/5.0))
        for j in range(nl):
            line=f.readline()
            a=line.strip().strip("(").strip(")").split(")  (") 
            b=[ (item.split('->')) for item in a]
            mp.extend(b)
        rotmap.append(mp)
        f.readline() 
    symop.append((rot, gtrans, rotmap, ptrans))
    f.readline() 

f.close()

# rot map of atoms in wannier projection
wann_atom_rotmap=np.zeros((nsymm,nptrans,natom), dtype=np.int32)
for isymm in range(nsymm):
    print "isymm", isymm, "---------"
    for iptr in range(nptrans): 
        print "iptrans", iptr, symop[isymm][3][iptr]
        for ii in range(natom): 
            target = np.int32(symop[isymm][2][iptr][ atom_num[ii]-1 ][1])
            indx=atom_num.index(target)
            wann_atom_rotmap[isymm,iptr,ii]=indx
            print "OUTCAR atom rot map:", atom_num[ii], "->", target, "=====>", "wannier atom rot map:", ii, "->", indx

if __name__=='__main__':
    f=open('symmetry.dat','w')
    for i in range(nsymm):
        f.write("rot: "+str(i+1)+"\n")
        rot=symop[i][0]
        for j in range(3):
            line="   ".join(rot[j])+"\n"
            f.write(line)
        gtrans=symop[i][1]
        line="   ".join(gtrans)+"\n"
        f.write(line)
        for iptr in range(nptrans):
            ptrans_tmp=symop[isymm][3][iptr]
            line="   ".join(ptrans_tmp)+"\n"
            f.write(line)
            rotmap=symop[i][2][iptr]
            for j in range(natom):
                line="   ".join(rotmap[j])+"\n"
                f.write(line)

            f.write("\n\n")
    f.close()


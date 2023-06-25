#!/bin/bash

# set global python path
export PYTHONPATH=/home/yuec/Source/wannier_tools/wannhr_symm_Mag/lib_Mag/:$PYTHONPATH

# step 1, perform symmetrization for magnetic operation without T
python2  /home/yuec/Source/wannier_tools/wannhr_symm_Mag/symmhr_addrptblock/symmhr_addrptblock_magnetic_conjugate_when_primed.py  1.0 wannier90_hr.dat
# after the first step, there will a file written out, which will be the input file for step 2. 

# step 2, perform symmetrization for magnetic operation with T
python2  /home/yuec/Source/wannier_tools/wannhr_symm_Mag/symmhr_addrptblock/symmhr_addrptblock_magnetic_conjugate_when_primed.py -1.0 wannier90_hr.dat_unprimed_nsymm_4



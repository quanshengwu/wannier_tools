Copyright: Changming Yue, Department of Physics, Southern University of Science and Technology, Shenzhen 518055, China. 
Email: yuecm@sustech.edu.cn; yuechangming8@gmail.com

# it is well tested with phonopy==1.11.8

A python tool based on Phononpy to generate phonon tight-binding Hamiltonian
whose format is similar to wannier90_hr.dat produced by Wannier90.

Author: Changming Yue, Institute of Physics, CAS.
email: yuechangming8@gmail.com

Necessary file: 
1. POSCAR: primitive cell atom positions used by vasp
2. FORCE_CONSTANTS: force constants matrix produced by vasp, abinit, etc.
3. input.in: to specify super cell dimenstion, reading force constants, etc
   Here is one sample:
             ATOM_NAME = Si Zr O
             DIM = 3 3 3
             FORCE_CONSTANTS = READ

Usage: python phonon_hr.py -d --dim="3 3 3" -c POSCAR -p band.conf
where dim="3 3 3" is just for a supercell with 3*3*3 size. You can set your
own size as you want.

Good luck!


About the installation:
1. We use Python2.7 which is not compatible with Python3
2. Please install phonopy before using this package.
3. You must have some knowledge of python 
4. Please add the current folder path to the system path and pythonpath like
export PATH=~/wannier_tools/phonopyTB:$PATH
export PYTHONPATH=~/wannier_tools/phonopyTB:$PYTHONPATH

5. If you don't set the path as in 4, please copy all the files into your working directory.

6. Please check the FeSi example

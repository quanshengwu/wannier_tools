# phonon_hr.py is well tested with phonopy-1.11.8

# In this folder we use WannierTools to perform the phonon calculation.
0. Quick running, to get the results in PRL
tar -xzvf phonopyTB_hr.dat.tar.gz
mpiexec -np 4 wt.x &

1. Firstly, you have to follow the instruction to install the phonopyTB in the root folder of WannierTools

2. Then prepare the necessary files to generate the TB model with the FORCE_CONSTANTS or FORCE_SETS.
A. POSCAR: primitive cell atom positions used by vasp
B. FORCE_CONSTANTS: force constants matrix produced by vasp, abinit, etc.
C. band.conf: to specify super cell dimenstion, reading force constants, etc
   Here is one sample:
   ATOM_NAME = Fe Si
   DIM = 3 3 3
   FORCE_CONSTANTS = READ
   BAND = 0.5 0 0 0 0 0 0.5  0.5  0.5  0.5  0.5 0 
   BAND_POINTS=100
   FC_SYMMETRY = 1
   fc_spg_symmetry = .true.

Usage: python phonon_hr.py --dim="3 3 3"  -p band.conf -c POSCAR
where dim="3 3 3" is just for a supercell with 3*3*3 size. You can set your
own size as you want.

Before using phonopyTB pacage, please be sure that you are familiar with the phonopy package. If you 
ecounter some problems with phonopyTB, most of them should be the settings for phonopy. 

Good luck!


#> in this folder
From FORCE_SETS to FORCE_CONSTANTS
phonopy --writefc -c POSCAR --dim='3 3 3'
phonon_hr.py -p band.conf -c POSCAR

Please add two lines below in the band.conf if consider the symmetry
FC_SYMMETRY = 1
fc_spg_symmetry = .true.

3. Most of the input tags in the wt.in are the same as the electron system. 

4. About the generation of Born charge and the dielectric tensor, please check the website of phonopy.



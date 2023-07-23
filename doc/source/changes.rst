
**Updates**
^^^^^^^^^^^^^^^^

.. role:: blue

**WannierTools V2.7.0** Jul 24 2023

1. Add a new function to calculate anomalous Nerest Effect (ANE).

2. Add a new function to calculate spin Hall conductivity (SHC).

3. Add a python code (utility/wannhr_symm_Mag) to symmetrize the magnetic tight-binding model hr.dat.

4. Add a post-processing tool (useful_scripts/post_sigma_OHE/post_sigma_OHE.py) to calculate the resistivity. The format of output files of conductivity tensor was changed.

5. Fixed a bug in QPI calculation.

6. Fixed several bugs related to atomic units

**WannierTools V2.6.0** Feb 24 2020

1. Add a new function to calculate the Landau levels including the Hofstadter butterfly and Wannier diagram.
 
2. Extend the HmnR format from dense format which is defined by Wannier90 to a sparse format which only store the 
   non-zero entries of HmnR.


**WannierTools V2.5.0** Dec 9 2019

1. Add a new function to calculate ordinary magnetoresistance of a nonmagnetic metal or semimetal with given magnetic field. See ref. Magnetoresistance from Fermi surface topology, ShengNan Zhang, QuanSheng Wu, Yi Liu, and Oleg V. Yazyev, Phys. Rev. B 99, 035142 (2019)

2. Add a new function to calculate unfolded bands from a supercell calculation. 


**WannierTools V2.4.2** July 9 2019

1. Write out POSCAR-slab to help user to perform DFT calculations. 

2. Add a new funtion to calculate the projected spin texture for the bulk system with vacuum or without inversion symmetry. 
   BulkSpintext_calc = T.
   Add a new CARD called SELECTED_ATOMS to control the projection. 

**WannierTools V2.4.0** Sep 1 2018

1. We have a big update in this version, since Changming Yue put his symmetrization code into WannierTools. So now we can symmetrize the Wannier functions based tight binding model. This funcionality is included in the wannhr_symm/ folder. 

2. Fixed a bug about reading Rcut. 

3. Distinguish two dual surface in the slab band calculation. 

4. Update an example about calculating the mirror Chern number of ZrTe.

**Develop branch** June 27 2018
We added the test version of phonon system. Welcome to git clone the develop branch and test it.

.. code:: console

   git clone https://github.com/quanshengwu/wannier_tools.git


**WannierTools V2.3.0**

1. Fixed a bug.
2. Add Translate_to_WS_calc in the CONTROL namelist. This works for BulkFS_plane_calc.

**WannierTools V2.2.9**

1. Fixed several bugs.

2. Added two new functionalities:

a. Calculate the energy levels at given k points in the KPOINTS_3D card. This is called the point mode.

b. Calculate anomalous Hall conductivity (AHC).

**WannierTools V2.2.6**

1. Discard the Miller indicies

2. Discard the third vector in the SURFACE card. The surface plane is specified only by two lattice vectors sitting on it.



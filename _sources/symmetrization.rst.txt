Symmetrization of wannier90_hr.dat **New**  usesful for ir2tb
=============================================================


Here is a brief introduction of the symmetrization functionalities. Basically, wannhr_symm is an independent package based on Python2.7 
written by Changming Yue (yuechangming8@gmail.com). Although this package is very useful, its requirement is also very restrict. 

1. The Wannier functions must be atomic like, which means there should be very weak hybridization between orbitals. It's not possible to deal with 
   the sp2, sp3 like wannier orbitals. 

2. At present, it can only take care of the s, p, d, t2g, eg orbitals. p orbitals means all three orbitals ordered with pz, px, py. 
   d orbitals are ordered as  “dz2”, “dxz”, “dyz”, “dx2-y2”, “dxy”. Such situation is automatically satisfied with in VASP+Wannier90

3. The file name of tight binding Hamiltonian should be 'wannier90_hr.dat'

4. This symmetrization process is crucial to run “ir2tb” to order to obtain the correct irreducible representations, since the Wannier90 usually does not enforce the symmetries of the system.


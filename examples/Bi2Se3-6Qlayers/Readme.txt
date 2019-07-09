example: Bi2Se3-6Qlayers

6 Quintuple-layer + 12 Angstrom vacuum Bi2Se3, We treated it as a bulk system. 
The spintexture of the surface states can be obtained by setting Bulkspintext_calc=T and
setting SELECTED_ATOMS such that the surface atoms are selected.

For example, in this case

SELECTED_ATOMS
2 ! number groups of selected atoms
6 12 18 24 30  ! top surface's atoms
1  7 13 19 25  ! bottom surface's atoms

There are two groups of selected atoms. One is for the top surface, the other one is for the bottom surface.
The indicies in the SELECTED_ATOMS are the same as the indicies in the ATOM_POSITIONS.

The outputs are bulkspintext.dat and bulkspintext.gnu. You can get the spintexture plot by run "gnuplot bulkspintext.gnu" 
for the first group. For other groups, you need to modify bulkspintext.gnu by choosing different columns that can be found
in the bulkspintext.dat.

The first 3 columns of bulkspintext.dat are the cartesian coordinates of k. The coordinate system is the same as LATTICE card in the wt.in.
The 4-6 columns are also the Cartesian coordinates of k, however, the coordinates system is a new one. 
The z axis of the new coordinate system is perpendicular to the kplane that defined in KPLANE_BULK. The x and y axis are sitting on that kplane.
The x axis is along the first vector in the kplane which is the second numerical line in the KPLANE_BULK.
The spin vectors are defined in the new coordinate system.

For bulk spin texture, you have to set the following necessary parameters:

Bulkspintext_calc=T 

Eta_Arc = 0.001     ! infinite small value, like brodening 
E_arc =  0.3         ! energy for calculate Fermi Arc
Nk1 = 101   ! number k points  odd number would be better
Nk2 = 101   ! number k points  odd number would be better

SELECTED_ATOMS
2 ! number groups of selected atoms
6 12 18 24 30  ! top surface's atoms
1  7 13 19 25  ! bottom surface's atoms


KPLANE_BULK
 0.00  0.00  0.00   ! center of 3D k plane 
 0.20  0.00  0.00   ! The first vector to define a k plane in 3D BZ
 0.00  0.20  0.00   ! The second vector to define a k plane in 3D BZ



# wannier_tools
Written by QuanSheng Wu in Fortran (wuquansheng@gmail.com)

Use wannier functions to get the surface state of slab system 
or edge state of nanowire system or just bulk band. Especially
usefull for topological insulator.

Jan 29 2015 at ETH
With given hr files from wannier90, 
1. we can get the bulk energy band with given 
high symmetry kpoint lines. 

2. We can get the energy band for slab system when 
you specify a plane, which can be done by modify ham_qlayer2qlayer.f90. From 
this calculation, you can get surface state. This is very useful for studying
topological insulator. 

3. We can get the energy band for nanowire system also.

4. We can get fermi-surface for slab system. 

5. We can use iterative green's 
function to get beautifull surface state.



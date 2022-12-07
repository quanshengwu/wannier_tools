# 2022.Dec.8
An example to calculate spin Hall conductivity of Pt
The DFT calculations are performed by using VASP with pseudopotential potpaw_PBE.54/Pt/POTCAR.
s,p,d orbitals of Pt are used as the projectors when construct Wannier functions using Wannier90 v1.2.

Necessary input parameters are included in the wt.in file. 

Run it.

tar xzvf wannier90_hr.dat.tar.gz
mpiexec -np 8 wt.x &

outputs: 27 components of sigma^{\gamma}_{\alpha, \beta}, where \alpha, \beta, \gamma= x, y, z
sigma_shc.txt and sigma_shc.gnu

# generate plots using gnuplot software
gnuplot sigma_shc.gnu

# It runs about 200s on Macbook pro m1max using 8 cores.

Benchmark:
sigma_xy^z(E=EF)=2166 (hbat/e)S/cm which is close to the one (2280 (hbat/e)S/cm) 
which is obtained from PHYSICAL REVIEW B 98, 214402 (2018)

Parameters used here:
Lattice constant: 3.9239 Angstrom
Kmesh : 41*41*41
Smearing : 0.05eV

# This is an example about twsited bilayer graphene (TBG)
# References: http://dx.doi.org/10.1103/PhysRevLett.126.056401
# https://pubs.acs.org/doi/10.1021/acs.nanolett.9b05117
# http://www.nature.com/articles/s41567-020-0825-9

# including obtaining band structure, Wilson loop and unfolded bands
# POSCAR-rigid: crystal structure of TBG with twist angle 7.34 degree.
# POSCAR-relaxed: LAMMPS relaxed crystal structure with CC.KC and CH.airebo potential

!> You can generate your hr.dat with the given tgtbgen program
# dense format storage of tight-binding parameters
# system.in-dense: an input file for tgtbgen program to generate the TG_hr.dat and wt.in
$ cp system.in-dense system.in
$ cp POSCAR-relaxed POSCAR
$ tgtbgen &

# dense format storage of tight-binding parameters
# system.in-dense: an input file for tgtbgen program to generate the TG_hr.dat and wt.in
$ cp system.in-dense system.in
$ cp POSCAR-relaxed POSCAR
$ tgtbgen &


!> Or you can use the generated ones.
#TG_hr.dat-sparse.tar.gz  : sparse-format stored 
#TG_hr.dat-dense.tar.gz   : dense-format stored 
tar xzvf TG_hr.dat-sparse.tar.gz 
tar xzvf TG_hr.dat-dense.tar.gz 



# obtain band structure with TG_hr.dat in dense format 
$ cp wt.in-bands-dense wt.in
$ mpirun -np 4 wt.x &
$ gnuplot bulkek.gnu
$ evince bulkek.pdf

# obtain band structure with TG_hr.dat in sparse format 
$ cp wt.in-bands-sparse wt.in
$ mpirun -np 4 wt.x &
$ gnuplot bulkek.gnu
$ evince bulkek.pdf

# obtain Wilson loop of the four 'flat' bands around the Fermi level
$ cp wt.in-wilsonloop-dense wt.in
$ mpirun -np 4 wt.x &
$ gnuplot wcc.gnu0
$ evince wcc.pdf

# obtain unfolded band along G-K-M-G path in the Graphene Brillouin zone
$ cp wt.in-unfolding-kpath wt.in
$ mpirun -np 4 wt.x &
$ gnuplot spectrum_unfold_kpath.gnu
$ eog spectrum_unfold_kpath.png

# obtain unfolded iso-band spectrum at an fixed energy E_arc around K point of the Graphene Brillouin zone
$ cp wt.in-unfolding-kplane wt.in
$ mpirun -np 4 wt.x &
$ gnuplot spectrum_unfold_kplane.gnu0
$ eog spectrum_unfold_kplane.png


# spectrum_unfold_kplane-kmesh-201x201.png: an unfolded spectrum with dense kmesh Nk1*Nk2

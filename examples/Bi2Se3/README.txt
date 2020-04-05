# Study bulk band structure, Z2 topological number and surface states of 3D strong topological insulator Bi2Se3

1.  Preparation: unzipp the tight-binding hamiltonian
    $ tar xzvf wannier90_hr.dat.tar.gz

2.  Run WannierTools:
    $ mpirun -np 2 wt.x &

3.  Plot bulk band structure
    $ gnuplot bulkek.gnu
    $ evince bulkek.pdf

4.  Check the topological number
    $ gnuplot wanniercenter3D_Z2.gnu-tutorial
    $ evince wanniercenter3D_Z2.pdf
    $ sed -n '/# z2 number/,/Time/p' WT.out

5.  Plot surface state spectrum and its spin-texture
    $ gnuplot surfdos_l.gnu
    $ gnuplot arc_l.gnu
    $ gnuplot spintext_l.gnu
    $ gnuplot surfdos_r.gnu
    $ gnuplot arc_r.gnu
    $ gnuplot spintext_r.gnu

6.  Plot band structure of slab system defined by SURFACE and Nslab
    $ gnuplot slabek.gnu
    $ eog slabek.png

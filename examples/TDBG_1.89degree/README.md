# AB-AB Twisted double bilayer graphene with twist angle 1.89 degree.

## Obtain tight-binding model. TDBG_m17_hr.dat.tar.gz

```
tar xzvf TDBG_m17_hr.dat.tar.gz
```

## Obtain the band structure

```
cp wt.in-bands wt.in
mpiexec -np 4 wt.x &
gnuplot bulkek.gnu
evince bulkek.pdf
```


## Obtain the Hofstadter-butterfly and Wannier diagram

```

cp wt.in-landaulevel wt.in
mpiexec -np 4 wt.x &
gnuplot LandauLevel_B_dos.gnu
gnuplot wannierdiagram.gnu
```

**Note** : Please increase Nslab gradually. The size of the Hamiltonian is equal to 
Nslab*Num_wann. The time of diagonalization the Hamiltonian is quadratically increase as Nslab.


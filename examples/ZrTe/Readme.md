# Obtain the Wannier tight-binding model

Please decompress ZrTe_hr.tar.gz to get wannier90_hr.dat
`tar xzvf ZrTe_hr.tar.gz`


# Berry curvature calculation
```
cp wt.in-berrycurvature wt.in

wt.x & 
or run it with MPI
mpiexec -np 4 wt.x &
```

Plot Berry curvature
`gnuplot Berrycurvature.gnu`


you will have Berrycurvature.png and Berrycurvature_EF.png


# Mirror Chern number calculation

```
cp wt.in-mirrorchernnumber wt.in
wt.x & 
or run it with MPI
mpiexec -np 4 wt.x &
```

a. Get mirror chern number from the WT.out by 
`grep 'MCN' WT.out`

b. plot the WCC/Wilson loop with 
`gnuplot wcc-mirrorchernnumber.gnu`


# Find node points (includes Weyl points and nodal lines)
```
cp wt.in-findnodes wt.in

wt.x & 
or run it with MPI
mpiexec -np 4 wt.x &
```

Information of all nodes is included in Nodes.dat.
All the nodes in the Nodes.dat are in the Wigner-Seitz cell of reciprocal lattice.

You can obtain a rough plot using gnuplot. You can use Pymatgen to generate the BZ information and plot them together.

`gnuplot Nodes.gnu`

# Now we are going to unfold the energy bands from the Brillouin zone of twisted bilayer graphene 
to the Brillouin zone of mono-layer Graphene

# Get more information at https://www.wanniertools.org/examples/band-unfolding/

1. Uncompress the tight binding Hamiltonian file
   $ tar xzvf tbg_hr.dat.tar.gz

2. Get band structure of twisted bilayer graphene
   $ cp wt.in-bands wt.in
   $ mpirun -np 2 wt.x 
   $ gnuplot bulkek.gnu
   $ evince bulkek.pdf

3. Get the unfolded bands
   $ cp wt.in-unfold_bands wt.in
   $ mpirun -np 2 wt.x 
   $ gnuplot spectrum_unfold_kpath.gnu
   $ eog spectrum_unfoldz_kpath.png

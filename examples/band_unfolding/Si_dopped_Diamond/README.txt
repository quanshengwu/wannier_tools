# Unfold the energy band of Si-dopped Diamond

1. Get band structure of Pure Diamond
   $ tar xzvf Diamond_hr.dat.tar.gz
   $ cp wt.in-Diamond wt.in
   $ wt.x
   $ gnuplot bulkek.gnu
   $ cp bulkek.pdf Diamond-bulkek.pdf


2. Get band structure of Si-dopped Diamond
   $ tar xzvf Si_dopped_Diamond_hr.dat.tar.gz
   $ cp wt.in-bands_Si_dopped_Diamond wt.in
   $ wt.x
   $ gnuplot bulkek.gnu
   $ cp bulkek.pdf Si_dopped_Diamond-bulkek.pdf

3. Unfold band structure of Si-dopped Diamond
   $ tar xzvf Si_dopped_Diamond_hr.dat.tar.gz
   $ cp wt.in-unfold_Si_dopped_Diamond wt.in
   $ wt.x
   $ gnuplot spectrum_unfold_kpath.gnu
   $ eog spectrum_unfoldz_kpath.png



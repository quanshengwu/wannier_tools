# Unfold the energy bands of C-vacancy Diamond

1. Get band structure of Pure Diamond
   $ tar xzvf Diamond_hr.dat.tar.gz
   $ cp wt.in-Diamond wt.in
   $ wt.x
   $ gnuplot bulkek.gnu
   $ cp bulkek.pdf Diamond-bulkek.pdf


2. Get band structure of C-vacancy Diamond
   $ tar xzvf C_vacancy_Diamond_hr.dat.tar.gz 
   $ cp wt.in-bands_C_vacancy_Diamond wt.in
   $ wt.x
   $ gnuplot bulkek.gnu
   $ cp bulkek.pdf C_vacancy_Diamond-bulkek.pdf

3. Unfold band structure of C-vacancy Diamond
   $ tar xzvf C_vacancy_Diamond_hr.dat.tar.gz
   $ cp wt.in-unfold_C_vacancy_Diamond wt.in
   $ wt.x
   $ gnuplot spectrum_unfold_kpath.gnu
   $ eog spectrum_unfoldz_kpath.png



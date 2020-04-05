# Here we start to study a 3D Weyl semimetal TaAs

1. Get tight binding model
   $ tar xzvf wannier90_hr.dat.tar.gz

2. calculate band structure
   $ cp wt.in-bands wt.in
   $ mpirun -np 2 wt.x &

3. Find all Weyl points
   $ cp wt.in-findnodes wt.in
   $ mpirun -np 2 wt.x &


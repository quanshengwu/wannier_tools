# This is an example to study Haldane model with three phases

1. Trivial insulator phase

1). Generate tight-binding model and input file for WannierTools
    $ python haldane_hr_gen-trivial-insulator.py
    $ cp wt.in-trivial-insulator  wt.in

2). Run WannierTools
    $ mpirun -np 2 wt.x &

3). Plot band structure
    $ gnuplot bulkek.gnu
    $ gnuplot bulkek_plane.gnu
    $ evince bulkek.pdf
    $ eog bulkek_plane.png

4). Plot Berry curvature
    $ gnuplot Berrycurvature.gnu 
    $ eog Berrycurvature.png

5). Plot Wannier charge center
    $ gnuplot wcc.gnu
    $ evince wcc.eps

6). Plot band structure of slab system
    $ gnuplot slabek.gnu
    $ eog slabek.png


2. Chern insulator phase

Do the same step as Trivial insulator phase but with haldane_hr_gen-chern-insulator.py
and wt.in-chern-insulator


3. Gapless semimetal phase

1). Run WannierTools 
    $ python haldane_hr_gen-gapless.py
    $ cp wt.in-gapless wt.in
    $ mpirun -np 2 wt.x &

2). Plot band structure
    $ gnuplot bulkek.gnu
    $ gnuplot bulkek_plane.gnu
    $ evince bulkek.pdf
    $ eog bulkek_plane.png

3). Plot band structure of slab system
    $ gnuplot slabek.gnu
    $ eog slabek.png



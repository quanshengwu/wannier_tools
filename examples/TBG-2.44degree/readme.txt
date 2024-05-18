# this is an example to calculate the 4 flat bands and their Wilson loop

cp wt.in-wilsonloop-flatbands wt.in
tar xzvf TBG_hr.dat.tar.gz

mpiexec -np 32 wt.x &

gnuplot wcc.gnu
gnuplot bulkek.gnu

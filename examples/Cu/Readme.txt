0. unzip the hr.dat.tar.gz
$ tar xzvf wannier90_hr.dat_nsymm48.tar.gz

1. Calculate band structure and Fermi surface
$ mkdir band
$ cp wt.in-bands band/wt.in
$ cp wannier90_hr.dat_nsymm48 band/
$ cd band
$ mpirun -np 4 wt.x &
$ gnuplot bulkek.gnu
$ xcrysden --bxsf FS3D.bxsf 

2. Calculate MR for different magnetic fields along Theta=0, 18, 30, 45 degree and Phi=90.
$ mkdir OHE-theta0
$ cp wannier90_hr.dat_nsymm48 OHE-theta0/
$ cp wt.in-OHE-theta0 OHE-theta0/wt.in
$ cd OHE-theta0
$ mpirun -np 4 wt.x &
$ cp rho_total_mu_0.00eV.dat ../rho_total_mu_0.00eV.dat-theta0
$ cd ..

repeat......

$ mkdir OHE-theta45
$ cp wannier90_hr.dat_nsymm48 OHE-theta45/
$ cp wt.in-OHE-theta0 OHE-theta45/wt.in
$ cd OHE-theta45
$ mpirun -np 4 wt.x &
$ cp rho_total_mu_0.00eV.dat ../rho_total_mu_0.00eV.dat-theta45
$ cd ..
$ gnuplot rhotheta.gnu

3. Get the bulk Fermi surface in a fixed k plane
$ mkdir FS-contour
$ cp wannier90_hr.dat_nsymm48 FS-contour/
$ cp wt.in-FS-contour FS-contour/wt.in
$ cd FS-contour
$ mpirun -np 4 wt.x &
$ gnuplot fs_kplane.gnu

4. Get the evolution of momentum k under magnetic field
$ mkdir OHE-evolve
$ cp wannier90_hr.dat_nsymm48 OHE-evolve/
$ cp wt.in-OHE-evolve OHE-evolve/wt.in
$ cd OHE-evolve
$ mpirun -np 4 wt.x &
Then open evolve_band_6.dat file and check whether there are k points evolution.

5. Calculate AMR
$ cd AMR-xy
$ cp ../wannier90_hr.dat_nsymm48 ./
$ chmod 777 xyplane.sh
replace the content between "cat>\$dir/wt-theta.sh<<EOF" and "EOF" by your own shell script and then run the xyplane.sh script:
$ ./xyplane.sh
when calculations are finished, run the python script and gnuplot script
$ python AMR_rhoT.py
$ cp AMR_rhotheta.py rho/
$ cd rho/
$ python AMR_rhotheta.py
$ cd rhotheta/
$ cp ../../AMR_rhoxx.gnu ../../AMR_rhozz.gnu ./
$ gnuplot AMR_rhoxx.gnu
$ gnuplot AMR_rhozz.gnu

6.you can set different relaxation times for different bands in post_sigma_OHE.py 

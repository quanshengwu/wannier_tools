# Added on Sep.05.2019 By QuanSheng Wu
This is a example for calculating ordinary magnetoresistance with given magnetic field direction.
1. Calculate band structure and Fermi surface
   $ cp wt.in-bands wt.in
   $ mpiexec -np 4 wt.x&
   $ gnuplot bulkek.gnu
   $ xcrysden --bxsf FS3D.bxsf &

2. Calculate MR for different magnetic fields along Theta=0, 18, 30, 45 degree and Phi=90.
   $ cp wt.in-OHE-theta0 wt.in
   $ mpiexec -np 4 wt.x &  
   $ cp  sigma_band_6_mu_0.00eV_T_30.00K.dat sigma_band_6_mu_0.00eV_T_30.00K.dat-theta0
   $ cp wt.in-OHE-theta18 wt.in
   $ mpiexec -np 4 wt.x &  
   $ cp  sigma_band_6_mu_0.00eV_T_30.00K.dat sigma_band_6_mu_0.00eV_T_30.00K.dat-theta18
   $ cp wt.in-OHE-theta30 wt.in
   $ mpiexec -np 4 wt.x &  
   $ cp  sigma_band_6_mu_0.00eV_T_30.00K.dat sigma_band_6_mu_0.00eV_T_30.00K.dat-theta30
   $ cp wt.in-OHE-theta45 wt.in
   $ mpiexec -np 4 wt.x &  
   $ cp  sigma_band_6_mu_0.00eV_T_30.00K.dat sigma_band_6_mu_0.00eV_T_30.00K.dat-theta45
   $ gnuplot rho.gnu
   $ evince rho.pdf


3. Get the evolution of momentum k under magnetic field
   $ cp wt.in-OHE-evolve wt.in
   $ mpiexec -np 4 wt.x &
   Then open evolve_band_6.dat file and check whether there are k points evolution.
   

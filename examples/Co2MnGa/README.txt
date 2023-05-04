# Here we start to calculate the temperature-dependent and chemical potential dependent anomalous Nernst coefficient of Co2MnGa
# Do not set Eta_arc too small, or there will be numerical instability.  

1. Get tight binding model
    $ tar xzvf wannier90_hr.dat.tar.gz

2. calculate anomalous Nernst coefficient with control tag ANE_calc = T. 
   The temperature ranges from TMin=20 to TMax=320 with NumT=33.
   The chemical potential ranges from E_fermi+OmegaMin to E_fermi+OmegaMax with OmegaNum=901.
   We can obtain the result written in ane.txt and the script for plotting in ane.gnu. 
   
    $ mpirun -np 2 wt.x &

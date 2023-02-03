# Here we start to calculate the temperature-dependent and chemical potential dependent anomalous Nernst coefficient of Co2MnGa
# With control tag ANEvsT = T or ANEvsEf = T, we can obtain the result written in ane.txt. 
# Do not set Eta_arc too small, or there will be numerical instability.  

1. Get tight binding model
    $ tar xzvf wannier90_hr.dat.tar.gz

2. calculate temperature-dependent anomalous Nernst coefficient with control tag ANEvsT = T. The temperature ranges from TMin to TMax.   
   Change E_FERMI if you want to change the chemical potential.  
   We give the reference result stored in ane.txt-varyT.  
   
    $ cp wt.in-varyT wt.in
    $ mpirun -np 2 wt.x &

3. calculate chemical potential dependent anomalous Nernst coefficient with control tag ANEvsT = F at 300 K by TMin=300, TMax=300, NumT=1.
   We calculate ANE ranging from (Ef-1, Ef+1) eV by OmegaMin=-1 and OmegaMax=1. 
   We give the reference result stored in ane.txt-varyEf.
    
    $ cp wt.in-varyEf wt.in
    $ mpirun -np 2 wt.x &

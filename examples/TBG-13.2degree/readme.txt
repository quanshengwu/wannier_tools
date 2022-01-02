This is an example to show how to unfold bands in WannierTools.
Because there are only 76 orbital in the unit cell, 
the tight-binding Hamiltonian is stored as a densed format 
which is the same format in Wannier90.


Please unzip the TBG_13.2_hr.tar.gz file in order to get TBG_13.2_hr.dat
tar xzvf TBG_13.2_hr.tar.gz

1. band calculation
cp wt.in-band wt.in
wt.x &

2. dos calcualtion 
cp wt.in-dos wt.in
wt.x &

3. band unfoldoing calcualtion 
cp wt.in-bandunfold wt.in
wt.x &


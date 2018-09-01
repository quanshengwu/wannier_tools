1. decompress ZrTe_hr.tar.gz to get wannier90_hr.dat
tar xzvf ZrTe_hr.tar.gz

2. run WannierTools as
wt.x &

3. plot Berry curvature
gnuplot  Berrycurvature.gnu

4. Get mirror chern number from the WT.out by 
grep 'MCN' WT.out

5. plot the WCC/Wilson loop with 
gnuplot wcc-mirrorchernnumber.gnu

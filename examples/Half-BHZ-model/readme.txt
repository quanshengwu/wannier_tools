# An example to show how does the band inversion and hybridization affect the band topology
# we prepared four python scripts and one wt.in file to calculate the band structure, slab bands
# AHC, Wilson loop, Berry curvature

# to run it, first, you need to generate the wannier TB model HalfBHZ_hr.dat
python3 Half_BHZ_hr_gen-case1.py

# 2nd run wanniertools wt.x
wt.x 

# 3rd get the plots

gnuplot bulkek.gnu
gnuplot slabek.gnu
gnuplot wcc.gnu
gnuplot sigma_ahc.gnu
gnuplot Berrycurvature.gnu

# use ll -tr to check the latest files
ll -tr


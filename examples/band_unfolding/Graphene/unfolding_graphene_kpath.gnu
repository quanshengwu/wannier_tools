#set terminal  postscript enhanced color font ",30"
#set output 'spectrum_unfold.eps'
set terminal pngcairo enhanced color font ",60" size 1920,1680
set palette defined ( 0 "white", 1  "#D72F01" )
set output 'spectrum_unfold_kpath.png'
set style data linespoints
set size 0.9, 1
set origin 0.05,0
unset key
set border lw 3
set view map
#set xtics font ",24"
#set ytics font ",24"
#set ylabel font ",24"
#set ylabel offset 1.5,0
emin=  -13.000000
emax=    7.000000
set xrange [0:    5.09086]
set ylabel "Energy (eV)"
set yrange [ emin : emax ]
set xtics ("M  "    0.00000,"K  "    0.84848,"G  "    2.54543,"K  "    4.24239,"M  "    5.09086)
set arrow from    0.84848, emin to    0.84848, emax nohead front lw 3
set arrow from    2.54543, emin to    2.54543, emax nohead front lw 3
set arrow from    4.24239, emin to    4.24239, emax nohead front lw 3
set colorbox
set cbrang[0:100]
set pm3d interpolate 2,2
splot 'spectrum_unfold_kpath.dat' u 1:2:(($8 *50))  w pm3d

set terminal pngcairo enhanced color font ",60" size 1920,1680
set palette defined ( 0 "white", 1  "#D72F01" )
set output 'spectrum_unfold_kplane.png'
set size 0.9, 1
set origin 0.05,0
set border lw 3
set pm3d
unset key
set view map
#set xtics font ",24"
#set ytics font ",24"
#set ylabel font ",24"
#set ylabel offset 1.5,0
set size ratio -1
set colorbox
set cbrange [0:50]
set pm3d interpolate 2,2
splot 'spectrum_unfold_kplane.dat' u 4:5:($13*50)  w pm3d

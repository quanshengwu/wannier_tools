#set terminal  postscript enhanced color font ",30"
#set output 'bulkek.eps'
set terminal pdf enhanced color font ",20"
set palette defined (-10 "#194eff", 0 "green", 10 "red" )
set output 'bulkek.pdf'
set style data points
set size 0.9, 1
set origin 0.05,0
unset key
set pointsize 0.1
#set xtics font ",24"
#set ytics font ",24"
#set ylabel font ",24"
#set ylabel offset 1.5,0
emin=   -2.385127
emax=    2.454274
set xrange [0:    0.66339]
set ylabel "Energy (eV)"
set yrange [ emin : emax ]
set xtics ("M  "    0.00000,"K  "    0.14019,"G  "    0.42057,"M  "    0.66339)
set arrow from    0.14019, emin to    0.14019, emax nohead
set arrow from    0.42057, emin to    0.42057, emax nohead
plot 'bulkek.dat' w p pt 7 ps 0.3 lc rgb 'black'

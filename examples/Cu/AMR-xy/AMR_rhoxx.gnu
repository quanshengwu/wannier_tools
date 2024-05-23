set encoding iso_8859_1
set terminal pdfcairo enhanced color font "Arial,20" size 8, 6
set output 'rhoxx.pdf'
set border lw 10
set autoscale fix
set ylabel '{/Symbol r}_{xx}*{/Symbol t} ({/Symbol W}*m*s)'
set format y "%1.1e"
set xlabel '{/Symbol \161}'
set xrange [0:180]
#set yrange [1.2e-21:2.5e-21]
set xtics 30
set key outside
#unset key
set palette defined (0 'red', 1 'green')
unset colorbox

set lmargin at screen 0.25    # 左边距占页面宽度的10%
set rmargin at screen 0.80    # 右边距占页面宽度的90%
set bmargin at screen 0.2    # 下边距占页面高度的10%
set tmargin at screen 0.9    # 上边距占页面高度的90%

set ylabel  offset -5,0
set xlabel  offset 0,-1
set key right vertical spacing 5
set xtics font ",33"
set ytics font ",33"
set ylabel font ",40"
set xlabel font ",40"
set key font ",27"
set key spacing 1.2
set key samplen 1.0
set xtics 30
set ytics 0.3e-21


Bmin =  2.00
Bmax = 10.00
NumB =   6
lw =    10

plot for [i=0:NumB-1] '90.0000K_Btau.dat' every :::i::i+1 u 1:2 w l lw lw  title sprintf('B=%.0f T',2*i)

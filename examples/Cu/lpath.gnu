set encoding iso_8859_1
set terminal pdfcairo enhanced color font "Times-New-Roman,30" size  8,8
set palette defined ( 0  "green", 5 "yellow", 10 "red" )
set output 'levolve.pdf'
set style data linespoints
set lmargin at screen 0.25    # 左边距占页面宽度的10%
set rmargin at screen 0.95    # 右边距占页面宽度的90%
set bmargin at screen 0.2    # 下边距占页面高度的10%
set tmargin at screen 0.9    # 上边距占页面高度的90%
unset ztics
set pointsize 0.8
set view 0,0
set key font ",40"
set key samplen 0.8
set key spacing 1.5
#set mxtics 0.005
#set ytics 0.05
set xtics font ",55"
set ytics font ",55"
set ylabel font ",60"
set xlabel font ",60"
set ylabel offset -3.8,0
set xtics  offset 0.0,-0.1
set xlabel offset 0.0,-0.5
set border lw 10


set xlabel "l_x ({\305})"
set ylabel "l_y ({\305})"
set xrange [-0.6:0.6]
set yrange [-0.6:0.6]
set xtics 0.3
set ytics 0.3
plot 'K.txt' u ( $5):($6)  w p pt 6 ps 0.2 lc rgb '#712A7D'   title sprintf('lz=0 ({\305})'),\

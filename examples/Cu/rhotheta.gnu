set encoding iso_8859_1
set terminal pngcairo enhanced color font ",30"  size 960, 800
set output 'rhotheta.png'
set border lw 2
set autoscale fix
set key left samplen 0.8
set ylabel "{/Symbol r}_{xx}*{/Symbol t} ({/Symbol W}*m*s)"
set xlabel "B{/Symbol t} (T.ps)"
set format y '%.0e'
plot 'rho_total_mu_0.00eV.dat-theta0'  every :::0::0 u 1:2 w l lw 3 title '{/Symbol q}=0', \
     'rho_total_mu_0.00eV.dat-theta18' every :::0::0 u 1:2 w l lw 3 title '{/Symbol q}=18', \
     'rho_total_mu_0.00eV.dat-theta30' every :::0::0 u 1:2 w l lw 3 title '{/Symbol q}=30', \
     'rho_total_mu_0.00eV.dat-theta45' every :::0::0 u 1:2 w l lw 3 title '{/Symbol q}=45'


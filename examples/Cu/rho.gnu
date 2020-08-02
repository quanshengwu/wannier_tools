set encoding iso_8859_1
set terminal pngcairo enhanced color font ",30"  size 960, 800
set output 'rho.png'
set border lw 2
set autoscale fix
set key left samplen 0.8
set ylabel "{/Symbol r}*{/Symbol t} ({/Symbol W}*m*s)"
set xlabel "B{/Symbol t} (T.ps)"
plot 'rho_band_6_mu_0.00eV_T_30.00K.dat-theta0' u 1:3 w l lw 3 title '{/Symbol q}=0', \
     'rho_band_6_mu_0.00eV_T_30.00K.dat-theta18' u 1:3 w l lw 3 title '{/Symbol q}=18', \
     'rho_band_6_mu_0.00eV_T_30.00K.dat-theta30' u 1:3 w l lw 3 title '{/Symbol q}=30', \
     'rho_band_6_mu_0.00eV_T_30.00K.dat-theta45' u 1:3 w l lw 3 title '{/Symbol q}=45'


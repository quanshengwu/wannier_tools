set encoding iso_8859_1
set terminal pngcairo enhanced color font ",30"  size 960, 800
set output 'rho.png'
set border lw 2
set autoscale fix
set key left samplen 0.8
set ylabel "{/Symbol r}*{/Symbol t} ({/Symbol W}*m*s)"
set xlabel "B{/Symbol t} (T.ps)"
set format y '%.1e'
plot 'rho_total_mu_0.00eV.dat-theta0' every :::0::0 u 1:2 w l lw 3 title '{/Symbol r}_{xx}*{/Symbol t}', \
     'rho_total_mu_0.00eV.dat-theta0' every :::0::0 u 1:6 w l lw 3 title '{/Symbol r}_{yy}*{/Symbol t}', \
     'rho_total_mu_0.00eV.dat-theta0' every :::0::0 u 1:10 w l lw 3 title '{/Symbol r}_{zz}*{/Symbol t}', \
     'rho_total_mu_0.00eV.dat-theta0' every :::0::0 u 1:3 w l lw 3 title '{/Symbol r}_{xy}*{/Symbol t}'


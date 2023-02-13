set terminal pdfcairo enhanced color font ",30" size 13, 6
set output 'ane.pdf'
set key outside
set palette defined (0 'red', 1 'green')
unset colorbox
##plot the temperature-dependent alpha_yx for the first six chemical potentials
#set xlabel "T (K)"
#set ylabel "{/Symbol a}_{yx} (A/(mK))"
#set ylabel offset 0.0,0
#plot for [i=0:5] 'ane.txt' every ::i*  33::i*  33+  32 u 2:(-$3) w l lt palette frac i/5. title sprintf('{/Symbol m}=%.3f eV',  -1.0+   2.0/ 900.0*i)

#plot the chemical potential dependent alpha_yx
set xlabel "T (K)"
set ylabel "{/Symbol a}_{yx} (A/(mK))"
set ylabel offset 0.0,0
plot for [i=0:  32] 'ane.txt' every   33::i u 1:(-$3) w l lt palette frac i/  32.0 title sprintf('T=%.3f K',  20.0+ 300.0/  32.0*i)

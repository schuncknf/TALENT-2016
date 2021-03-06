set xlabel "r(fm)"
set ylabel "density(fm-3)"
set style line 1 lt 2 lc rgb "red" lw 3
set linestyle 2 lt 5 lw 3
m="./density_n_new.dat"
n="./density_p_new.dat"
set terminal x11 0
set grid
set title 'Densities'
plot m using 1:2 w l ls 2 title 'Neutrons' 
replot n using 1:2 w l ls 1  title 'Protons'
set xrange[0:10]
replot

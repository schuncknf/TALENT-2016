set xlabel "r(fm)"
set ylabel "density(fm-3)"
m="./density_n_new.dat"
n="./density_p_new.dat"
set terminal x11 0
set nokey
set grid
set title 'Densities'
plot m using 1:2 with linespoints  
replot n using 1:2 with linespoints  

f = sprintf('occ_%f.dat',t)
set style line 1 lc rgb '#0060ad' lt 1 lw 0.6 pt 7 ps 0.8
set style fill transparent solid 0.5 noborder
set key off
set grid back linestyle 81
set xtics nomirror 
set ytics nomirror
set xrange[1.518:60]
set yrange[0:1]
set xlabel "Single-Particle Energies (MeV)"
set ylabel "Occupation probabilities"
plot f u 2:3 w filledcurves x1 ls 1,f u 2:3 w lp ls 1
t = t+0.08
pause 0.5
reread 

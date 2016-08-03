f = sprintf('density_%f.dat',t)
set style line 1 lc rgb '#0060ad' lt 1 lw 0.6 pt 7 ps 0.8
set style fill transparent solid 0.5 noborder
set key off
set grid back linestyle 81
set xtics nomirror 
set ytics nomirror
set xrange[0:6]
set yrange[0:1]
set xlabel "r (fm)"
set ylabel "Density (fm-3)"
plot f u 1:2 w filledcurves x1 ls 1,'density_hf' u 1:2 w l
t = t+0.04
pause 0.1
reread 

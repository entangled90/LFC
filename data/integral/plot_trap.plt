set border linewidth 1.5
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5;
set grid
set title 'Errore deterministici'
set xlabel 'log N'
set ylabel 'Errore'
unset key
set logscale x
set logscale y
set term postscript enhanced color landscape lw 1 'Verdana,10'
set output 'plot_trap.eps'
plot 'trap.dat' using 1:2 with linespoints, 'trap.dat' using 1:3 with linespoints
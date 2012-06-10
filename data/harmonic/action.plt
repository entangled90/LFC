reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
set xlabel 'N sweeps'
set ylabel 'S euclidea'
set term postscript enhanced color landscape lw 1 'Verdana,10'
set output 'action.eps'
set logscale x
plot 'action.dat' with linespoints

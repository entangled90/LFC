reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
set xlabel '|l-k|'
set ylabel 'Correlatore'
set term postscript enhanced color landscape lw 1 'Verdana,10'
set output 'harmonic.eps'
plot 'harmonic.dat' with linespoints

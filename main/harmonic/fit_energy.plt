reset
set border linewidth 1.5
set style line 1 lc rgb '#bf0d23' lt 1 lw 2 pt 7 ps 0.5 # --- red
set grid
set title "Metropolis"
set tics out nomirror
set xlabel "{/Symbol D}E"
set ylabel "Frequency"
set term postscript enhanced color landscape lw 1 "Verdana,10"
set output '../../data/harmonic/fit-energy.eps'
m= 0.962
s= 1e-3
pi=3.14159265
A = 1.5
n=100
f(x)=A*exp(-0.5*((x-m)/s)**2)/(2.50662827*s)
fit f(x) 'bin-energy.dat' using 1:2:3 via A,s,m
n=100#number of intervals
max=0.98;
min=0.943;
width=(max-min)/n#interval width
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
set xtics min,(max-min)/5,max
set boxwidth width
set style fill solid 0.5#fillstyle
set tics out nomirror
ti = sprintf("Gaussian Fit: \n { /Symbol m } = %f; {/Symbol s} = %f", m, s)
plot '../../data/harmonic/energy.dat' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb '#00ff00' title 'Binned data',f(x) w l ls 1 title ti
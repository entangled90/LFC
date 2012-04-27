#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>
/*
 *
 *
 *
 * INCOMPLETO!!!
 */
void plot(const char* data1,const char* data2,const char* data3 ){
    FILE *pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "reset\n");
    fprintf(pipe, "set border linewidth 1.5\n");
    fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n");
    fprintf(pipe, "set grid\n");
	fprintf(pipe, "set title 'Load Effect'\n");
    fprintf(pipe, "set xlabel 'N)'\n");
    fprintf(pipe, "set ylabel 'Valore integrale'\n");
    //fprintf(pipe, "unset key\n");
    fprintf(pipe, "set term postscript enhanced color landscape lw 1 'Verdana,10'\n");
    fprintf(pipe, "set output 'plot.eps'\n");
	fprintf(pipe, "plot '%s' title \"flat\" using 1:2:3 with errorbars w l ls 1,'%s' title \"gauss\" using 1:2:3 with errorbars w l ls 1,'%s' title \"root\" using 1:2:3 with errorbars w l ls 1 \n",data1,data2,data3);
    pclose(pipe);
}

#endif


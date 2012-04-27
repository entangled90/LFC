#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>
#include <stdlib.h>

void plot(const char* data1,const char* data2,const char* data3 ){
    FILE *pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "reset\n");
    fprintf(pipe, "set border linewidth 1.5\n");
    fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n");
    fprintf(pipe, "set grid\n");
	fprintf(pipe, "set title 'Campionamento importanza'\n");
    fprintf(pipe, "set xlabel 'N)'\n");
    fprintf(pipe, "set ylabel 'Valore integrale'\n");
   
    //fprintf(pipe, "unset key\n");
	fprintf(pipe,"set logscale x \n");
    fprintf(pipe, "set term postscript enhanced color landscape lw 1 'Verdana,10'\n");
    fprintf(pipe, "set output '../data/plot.eps'\n");
	fprintf(pipe, "plot '%s' with errorbars ,'%s' with errorbars  ,'%s' with errorbars  \n",data1,data2,data3);
	pclose(pipe);
}


void plot_noise(const char* data1,const char* data2,const char* data3 ){
    FILE *pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "reset\n");
    fprintf(pipe, "set border linewidth 1.5\n");
    fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n");
    fprintf(pipe, "set grid\n");
	fprintf(pipe, "set title 'Campionamento importanza'\n");
    fprintf(pipe, "set xlabel 'N'\n");
    fprintf(pipe, "set ylabel 'Rumore * sqrt(n)'\n");
    //fprintf(pipe, "unset key\n");
	fprintf(pipe,"set logscale x \n");
    fprintf(pipe, "set term postscript enhanced color landscape lw 1 'Verdana,10'\n");
    fprintf(pipe, "set output '../data/plot_noise.eps'\n");
	fprintf(pipe, "plot '%s' ,'%s'  ,'%s' \n",data1,data2,data3);
	pclose(pipe);
	}

void plot_harmonic (const char* data1, const char* out){
    FILE *pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "reset\n");
    fprintf(pipe, "set border linewidth 1.5\n");
    fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n");
    fprintf(pipe, "set grid\n");
	fprintf(pipe, "set title 'Harmonic'\n");
    fprintf(pipe, "set xlabel 'N'\n");
    fprintf(pipe, "set ylabel 'S'\n");
    fprintf(pipe, "unset key\n");
	/*
	 * fprintf(pipe,"set logscale x \n");
		fprintf(pipe,"set xrange[0:31] \n");
		fprintf(pipe,"set yrange[0:100] \n");
	*/
	fprintf(pipe, "set term postscript enhanced color landscape lw 1 'Verdana,10'\n");
	fprintf(pipe, "set output '%s'\n",out);
	fprintf(pipe, "plot '%s' with linespoints\n",data1);
	pclose(pipe);
}


#endif


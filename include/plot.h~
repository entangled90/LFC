#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "varie.h"

void plot(const char* data1,const char* data2,const char* data3 ){
    FILE *pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "reset\n");
    fprintf(pipe, "set border linewidth 1.5\n");
    fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n");
    fprintf(pipe, "set grid\n");
	//fprintf(pipe, "set title 'Campionamento importanza'\n");
    fprintf(pipe, "set xlabel 'N'\n");
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



void fit( const char * input ,const char * output, int cases ){
    double width;
    int* freq;
    int n = 100;
    int i ;
 //   double mean;
  //  double sigma;
    float tmp;
    double max;
    double min ;
    int n_events = 1;
    freq =  malloc(n*sizeof(int));
  /* EX funzione binning */
    for(i=0;i<n;i++)
      freq[i]=0;
    FILE* f=fopen(input,"r");
	fscanf(f,"%e\n",&tmp);
	max = tmp;
	min = tmp;
    while(fscanf(f,"%e\n",&tmp) == 1){
      n_events++;
		if(tmp>max)
			max = tmp;
		if (tmp < min)
			min = tmp;
    }
    fclose(f);
    f=fopen(input,"r");
    width=(max-min)/(double)n;
    while( fscanf(f,"%e\n",&tmp) == 1){
      for(i=0;i<n;i++){
	  if( (tmp>min+i*width) && (tmp<=min+(i+1)*width) ){
	  freq[i]++;
	  }
	  }
	  
    }
	fclose(f);
   switch(cases){
     case 1:
     f=fopen("bin-energy.dat","w");
     break;
     case 2:
     f=fopen("bin-matrix.dat","w");
   }
    for(i=0;i<n;i++){
        fprintf(f,"%lf\t%d\n",min+(i+0.5)*width,freq[i] );//,sqrt((double)freq[i]));
    }
    fclose(f);
    free(freq);
}


#endif


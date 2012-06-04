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



void fit( const char * input ,const char * output, double mean , double sigma){
    double width;
    int* freq;
	int n = 100;
	int i ;
	double tmp;
	double max;
	double min ;
	int n_events = 1;
	freq =  malloc(n*sizeof(int));
	/* EX funzione binning */
	for(i=0;i<n;i++)
        freq[i]=0;
    FILE* f=fopen(input,"r");
	fscanf(f,"%lf \n",&tmp);
	max = tmp;
	min = tmp;
    while(fscanf(f,"%lf \n",&tmp) == 1){
		n_events++;
		if(tmp>max)
			max = tmp;
		if ( tmp < min)
			min = tmp;
        
	}
	width=(max-min)/n,tmp;
	while(fscanf(f,"%lf \n",&tmp) == 1){
		for(i=0;i<n;i++)
			if(tmp>min+i*width && tmp<=min+(i+1)*width)
                freq[i]++;
	}
	fclose(f);
    f=fopen("bin.dat","w");
    for(i=0;i<n;i++)
		//if(freq[i]>4)
            fprintf(f,"%lf\t%d\n",min+(i+0.5)*width,freq[i]);
    fclose(f);
    free(freq);
	/* Fine ex funzione binning */
    FILE *pipe = popen("gnuplot -persist","w");
	fprintf(pipe, "reset\n");
	fprintf(pipe, "set border linewidth 1.5\n");
	fprintf(pipe, "set style line 1 lc rgb '#bf0d23' lt 1 lw 2 pt 7 ps 0.5 # --- red\n");
	fprintf(pipe, "set grid\n");
	fprintf(pipe, "set title \"Metropolis\"\n");
	fprintf(pipe, "set tics out nomirror\n");
	fprintf(pipe, "set xlabel \"{/Symbol D}E\"\n");
	fprintf(pipe, "set ylabel \"Frequency\"\n");
	fprintf(pipe, "set term postscript enhanced color landscape lw 1 \"Verdana,10\"\n");
	fprintf(pipe, "set output '%s'\n",output);
	fprintf(pipe, "m= %lf\n",mean);
	fprintf(pipe, "s= %lf\n", sigma);
	fprintf(pipe, "pi=3.14159265\n");
	//fprintf(pipe, "A = %lf \n", n_events*width);
	fprintf(pipe, "f(x)=exp(-0.5*((x-m)/s)**2)/(2.50662827*s)\n");
	//printf("A vale %lf \n ", n_events*width);
	fprintf(pipe, "fit f(x) 'bin.dat' via s,m\n");
	fprintf(pipe, "n=100\t#number of intervals\n");
	fprintf(pipe, "max= %lf \t#max value\n", max);
	fprintf(pipe, "min= %lf \t#min value\n",min);
	fprintf(pipe, "width=(max-min)/n\t#interval width\n");
	fprintf(pipe, "hist(x,width)=width*floor(x/width)+width/2.0\n");
	fprintf(pipe, "set xrange [min:max]\n");
	fprintf(pipe, "set yrange [0:]\n");
	fprintf(pipe, "set xtics min,(max-min)/5,max\n");
	fprintf(pipe, "set boxwidth width\n");
	fprintf(pipe, "set style fill solid 0.5\t#fillstyle\n");
	fprintf(pipe, "set tics out nomirror\n");
	fprintf(pipe, "ti = sprintf(\"Gaussian Fit:\\n { /Symbol m } = %%f; {/Symbol s} = %%f\", m, s)\n");
	fprintf(pipe, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb '#00ff00' title 'Binned data',f(x) w l ls 1 title ti\n",input);
	fclose(pipe);
    system("rm fit.log");
    system("rm bin.dat");
}


#endif


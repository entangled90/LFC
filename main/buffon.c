#define MAIN_PROGRAM
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "varie.h"
#include <time.h>
#include <stdlib.h>

#define L 0.99
#define PI 4*atan(1) 
/*
 * Global Variables
 */
	int N =1e7;
	int toss = 1e2;
	

double singlePi (double * x, double * y){
	double integral = 0;
	ranlxd(y,N);
	ranlxd(	x, N);
	int j =0;
	for (j = 0; j < N; j++){
		if (  (y[j] <  L*sin(PI*x[j])/2.0) || (y[j] > 1.0-L*sin(PI*x[j])/2.0) ){
			integral += 1;
		}
	}
	integral /= N;
	 return ( 2* L /integral);
}

long double meanPi (){
	
	int i = 0;
	double sum = 0;
	double *y;
	double *x;
	x = malloc(N * sizeof(double));
	y=malloc (N*sizeof(double));
	for (i = 0 ; i < toss ; i++){
		sum += singlePi(x,y)/toss;
		printf(" Siamo a  %d su %d \n", i+1, toss);
	}
	free(x);
	free(y);
	return( sum );
}

int main (int argc, char *argv[]) {
	double Pi  = 0;
	rlxd_init(1,time(NULL)%5000);
	Pi= meanPi();
	printf("Pi = %lf \n", (double) Pi);
	return (0);
	}

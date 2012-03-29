#define MAIN_PROGRAM
#include <stdio.h>
#include <math.h>
#include "random.h"
#include <stdbool.h>
#include "varie.h"
#include <time.h>
#include <stdlib.h>

#define L 0.99
#define NUM_OF(x) (sizeof (x) / sizeof *(x))
#define PI 4*atan(1) 
/*
 * Global Variables
 */
	static long long int N =1e7;
	static long int  toss = 1e1;
	

double singlePi (){
	double integral = 0;
	double *y;
	double *x;
	x = malloc(N * sizeof(double));
	y=malloc (N*sizeof(double));
	ranlxd(y,N);
	ranlxd(	x, N);
	 long int j =0;
	for (j = 0; j < N; j++){
		if (  (y[j] <  L*sin(PI*x[j])/2.0) || (y[j] > 1.0-L*sin(PI*x[j])/2.0) ){
			integral += 1;
		}
	}
	integral /= N;
	free(x);
	free(y);
	 return ( 2* L /integral);
}

long double meanPi (){
	
	long int i = 0;
	double sum = 0;
	for (i = 0 ; i < (long int) toss ; i++){
	sum += singlePi()/toss;
	if ( i % 10 == 0)
		printf(" Siamo a  %ld su %ld \n", i, toss);
	}
	return( sum );
}

int main (int argc, char *argv[]) {

	double Pi  = 0;
	Pi= meanPi();
	printf("Pi = %lf \n", (double) Pi);
	return (0);
	}

#define MAIN_PROGRAM
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "varie.h"
#include <time.h>
#include <stdlib.h>


#define N_NUMBERS 1e8
/* Lunghezza vettore nunmeri estratti*/
#define N  1e2

#define N_CYCLE  ( N_NUMBERS/N)


/*
 * Global Variables
 */
	double *y;
	double *x;

double singlePi (double * x, double * y){
	double integral = 0;
	ranlxd(y,N);
	ranlxd(	x, N);
	int j =0;
	for (j = 0; j < N; j++){
		if (  x[j]*x[j] + y[j]*y[j] < 1 ){
			integral += 1;
		}
	}
	integral /= N;
	 return ( 4*integral);
}

long double meanPi (){
	int i = 0;
	double sum = 0;
	x = malloc(N*sizeof(double));
	y = malloc (N*sizeof(double));
	for (i = 0 ; i < N_CYCLE ; i++){
	  sum += singlePi(x,y)/N_CYCLE;
	  if (  (int) ((i*100.0)/N_CYCLE)  != (int) ((i-1)*100.0/N_CYCLE) ) 
	    if( ((int) ((i*100.0)/N_CYCLE))% 10 == 0)
	      printf("%d % \n", (int)( i*100 / N_CYCLE)) ;
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

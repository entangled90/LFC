#define INTEGRATION_C
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"
#include "varie.h"

double Simpson ( double min , double max, double (*f) (double) ){
		return ( (max-min)*( f(min ) + 4.0*f(min + (max-min)/2 ) + f(max) )/6.0);
	}
	
double trapezio ( double min , double max, double (*f) (double) ){
	return ((max-min)*( f( min ) + f(max) )/2.0 ) ;
	}

double gaussianQuad ( double min , double max , double (*f) (double)){
	/*Viene usato un polinomio di grado 5 */
	double zero[5] = { 	0 ,
						sqrt(245.0 - 14.0*sqrt(70.0))/21.0,
						-sqrt(245.0 - 14.0*sqrt(70.0))/21.0,
						sqrt(245.0 + 14.0*sqrt(70.0))/21.0,
						-sqrt(245.0 + 14.0*sqrt(70.0))/21.0};
	double weight[5] = { 			(double)128.0/225.0,
						(double)1/900.0*(322.0 + 13*sqrt(70.0)),
						(double)1/900.0*(322.0 + 13*sqrt(70.0)),
						(double)1/900.0*(322.0 - 13*sqrt(70.0)),
						(double)1/900.0*(322.0 - 13*sqrt(70.0))} ;
	double integral =  0;
	/* Porto [min,max] in [-1,1] */
	int i = 0;
	for (i = 0 ; i< 5 ; i++){
		integral += weight[i]*f((max-min)/2.0*zero[i]+(max+min)/2.0);
	}
	return integral*(max-min)/2.0 ;
}


/*
 *  MethodNumber:
 * 	1 -> Trapezi
 *  2 -> Simpson
 *  3 -> Quadrature gaussiane
 * 	*f = integrand Function
 *  n = Number of sub intervals
 */

double partition (double min, double max , int n, int methodNumber , double (*f) (double) ){
	double h  = 0;
	int i = 0;
	double Sum = 0;
	//int swap = 0;
	double (*method) (double ,double , double (*) (double));
	switch(methodNumber){
		case 1:
			method = trapezio;
			break;
		case 2 :
			method = Simpson;
			break;
		case 3 :
			method = gaussianQuad;
			break;
		default:
			printf("Bad integration method!");
			exit(EXIT_FAILURE);
			break;
	}
	/*
	if ( max < min){
		h = min;
		min = max;
		max = h;
		swap = 1;
	}
	* */
	h = (max-min)/n;
	for (i = 0; i< n; i++){
		Sum += method( min + i*h, min + (1+i)*h, f);
		}
	return Sum;
	/*	
	if ( swap == 0)	
		return Sum;
	
	else if ( swap == 1) {
		return -Sum;
	}
	else{
		exit(1);
	}
	*/
}




#define VARIE_C
#include <math.h>
#include "varie.h"
#include "global.h"
static const double PI = 4*atan(1);

float meanOfFloatArray( float array[], int n){
	int i = 0;
	float sum = 0;
	for (i= 0; i< n; i++){
		sum += array[i];
		
		}
	return ( sum/n);
	}


	double meanOfDoubleArray( double array[] , int n){
	int i = 0;
	double sum = 0;
	for (i= 0; i< n; i++){
		sum += array[i];
		
		}
	return ( sum/n);
	}

	double flatPdf(double min, double max,double x){
		return  1/(max-min);
	}

	double gaussPdf(double variance, double mean, double x){
		return  exp(-(mean-x)*(mean-x)/(2*variance))/sqrt(2*PI*variance) ;
	}

	double gaussMomentPdf(double variance, double mean , int * n , double x){
	return pow(x,*n)/sqrt(2*PI*variance)*exp(-(mean-x)*(mean-x)/(2*variance));
	}
	
	double root_exp_pdf(double x){
	    return 2/sqrt(PI)*exp(-x)*sqrt(x);
	}

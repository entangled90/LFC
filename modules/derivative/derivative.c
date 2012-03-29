#define DERIVATIVE_C
#include <stdlib.h>
#include "utils.h"
#include <math.h>

double Derive( double (*function) (double ), double x, int n){
	double tmp = 0;
	double h = 0.1;
	if ( n< 0)
		exit(1)	; // ERROR!
	
	else if ( n == 0){	
		return( function(x)) ;
		}	
	else if ( n > 0){
		do{
			tmp= (Derive(function, x+h/2,n-1) - Derive(function,x-h/2,n-1))/h;
			h *= 0.5;
		} while ( fabs((Derive(function, x+h/2,n-1) - Derive(function,x-h/2,n-1))/h - tmp )> 1e-8 );
		return tmp;
	}
}

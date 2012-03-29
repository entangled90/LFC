#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "integration.h"
#include "varie.h"
#define MAIN_PROGRAM 

// #define N_SUB_INTERVAL 1000
#define integrandFunction2 log(1+x)
/*
 * asd 
*/
double parabola (double x ){
	return pow(x,2);
	}


double integrandLog( double x ){
	return ( log(1+x));
	}
double integrandPoly (double x ){
	return (pow(x,9) - pow(x,7) + 3.0);
	}
double dabs(double x){
	if (x >= 0)
		return x;	
	else 
		return (-x);
}

double D2Poly (x){
	return (72*pow(x,7) - 42*pow(x,5));	
	}





double primitiveLog (double x){
	return	(-x + log(1 + x) + x*log(1 + x)) ;	
	}
double primitivePoly (double x ){
		return	(pow(x,10)/10 - pow(x,8)/8 + 3*x) ;	
	}

double realIntegration (double (*primitive) (double x),double xmin, double xmax){
	return (primitive(xmax) - primitive(xmin));	
	}


int main (int argc, char *argv[]) {
	double xMin =  atof(argv[1]);
	double xMax =  atof(argv[2]);
	// double x = atof(argv[4]);
	int n =  atoi(argv[3]);
	int choice = 0;	
	double (* integrand ) (double  );
	double (* primitive) (double );
	
	printf("Scegli la funzione da integrare : (1) = Log , (2) = Poly \n");
	scanf("%d", &choice);
	switch(choice) {
		case 1 :
			integrand = (integrandLog);	 //, double xMin, double xMax);
			primitive = primitiveLog;
		break;
		case 2 : 
			integrand = (integrandPoly);
			primitive = primitivePoly;
		break;		
		default:
			printf("Inserisci 1 o 2 ");
			return (EXIT_FAILURE);
			break;
	}
	double integral_trap = partition(xMin,xMax, n, 1, integrand );
	double integral_simp = partition(xMin,xMax, n, 2, integrand );
	double integral_quad = gaussianQuad(xMin,xMax, integrand);
	double integral_true = realIntegration(primitive , xMin,xMax);
	printf("Valore esatto \t Trapezi \t Simpson \t gaussQuad\t Err. Trap \t Err. Simp \n");
	printf("%10lf \t %10lf \t %10lf \t %lf ", integral_true , integral_trap, integral_simp,integral_quad );
	printf("%e \t %e \n ", (integral_true - integral_trap) , (integral_true - integral_simp) );
	printf("\n");
	return(EXIT_SUCCESS);
}

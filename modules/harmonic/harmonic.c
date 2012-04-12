#ifndef HARMONIC_C
#define HARMONIC_C
#include <math.h>
#include <stdlib.h>
#include "harmonic.h"
#include "globalharmonic.h"
#include "random.h"

double potential (double x ){
	return (M*x*x*W*W/2.0);
	}

double lagrangian (double x2, double x1){
	return ( M/2.0*(x2-x1)*(x2-x1) -potential(x1)/2.0 - potential(x2)/2.0);
	}


double elagrangian (double x2, double x1){
	return ( M/2.0*(x2-x1)*(x2-x1) + potential(x1)/2.0 + potential(x2)/2.0);
	}


double efirst_action (double x[], int n){
	double S = 0.0;
	int i = 0;
	for(i=0; i<n ; i++){
		if(  i == n-1)
			S+= elagrangian(x[0],x[i]);
		else
			S += elagrangian(x[i+1],x[i]);
		}
	return S;
	}
double edelta_action ( double x_old[], double x_new, int position){
	double dS = 0.0;
	dS = elagrangian(x_old[(position+1)%Nx ], x_new) + elagrangian( x_new,x_old[(position-1+Nx)%Nx ])
		- elagrangian(x_old[(position+1)%Nx],x_old[position]) - elagrangian(x_old[position],x_old[(position-1+Nx)%Nx]);
	return dS;
	}


double 	first_action (double x[], int n){
	double S = 0.0;
	int i = 0;
	for(i=0; i<n ; i++){
		if(  i == n-1)
			S+= lagrangian(x[0],x[i]);
		else
			S += lagrangian(x[i+1],x[i]);
		}
	return S;
	}
double delta_action ( double x_old[], double x_new, int position){
	double dS = 0.0;
	dS = lagrangian(x_old[position+1], x_new) + lagrangian( x_new,x_old[position-1])
		- lagrangian(x_old[position+1],x_old[position]) - lagrangian(x_old[position],x_old[position-1]);
	return dS;
	}

void metropolis( double x[] , int position, double x_new){
	double dS= edelta_action(x,x_new,position);
	double tmp ;
	ranlxd(&tmp,1);
	if( tmp < exp(-dS) ){
		x[position] = x_new;
	}
}

double correlation ( double x[], double dK){
	int i = 0;
	double sum = 0.0;
	for(i = 0; i< Nx; i++){
		sum += x[i]*x[(i+dK)%Nx];
	}
	sum /= Nx
	return sum;
	}

#endif

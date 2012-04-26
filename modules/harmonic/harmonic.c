#ifndef HARMONIC_C
#define HARMONIC_C
#include <math.h>
#include <stdlib.h>
#include "harmonic.h"
#include "globalharmonic.h"
#include "random.h"
#include "struct.h"
#include <stdio.h>

double potential (double x ){
	return (M*x*x*W*W/2.0);
	}

double elagrangian (double x2, double x1){
	return ( M/2.0*(x2-x1)*(x2-x1) + potential(x1)/2.0 + potential(x2)/2.0);
	}


double efirst_action (double *x, int n){
	double S = 0.0;
	int i = 0;
	for(i=0; i<n ; i++){
			S += elagrangian(x[(i+1)%Nx],x[i]);
		}
	return S;
	}
double edelta_action ( double *x_old, double x_new, int position){
	double dS = 0.0;
	dS = elagrangian(x_old[(position+1)%Nx], x_new) + elagrangian( x_new,x_old[((position-1+Nx)%Nx)])
		- elagrangian(x_old[((position+1)%Nx)],x_old[position]) - elagrangian(x_old[position],x_old[((position-1+Nx)%Nx)]);
	return dS;
	}

void metropolis( double *x , int position, double *x_new){
	double dS= edelta_action(x,*x_new,position);
	double tmp ;
	ranlxd(&tmp,1);
	/* (1-dS + dS*dS/2 - dS*dS*dS/6) */
	if( tmp <  exp(-dS)){
		x[position] = *x_new;
	}
}

double  correlation ( double *x, int dK){
	int i = 0;
	double sum = 0.0;
	for(i = 0; i< Nx; i++){
		sum += x[i]*x[((i+dK)%Nx)];
	}
	sum /= Nx;
	return sum;
	}

	
	/*
	* Energy è un array di dimensione K_MAX - K_START.
	* input è un array di dimensione K_MAX
	*/
	
void DeltaE_cluster ( cluster_jk * input , cluster_jk *energy){
	int k,i;
	for( k = 0 ; k < K_MAX -K_START ; k++){
		for(i = 0; i< input->n_conf ; i++){
		energy[k].a[i] = acosh((input[(k)].a[i] + input[(k+2)].a[i])/( 2*input[k+1].a[i] ));
		}
	energy[k].mean  = acosh((input[(k)].mean + input[(k+2)].mean)/( 2*input[k+1].mean ));
	}
}

void matrix_element_cluster ( cluster_jk * input_E, cluster_jk * input_corr , cluster_jk *matrix_cluster){
	int k,i;
	for( k = 0 ; k < K_MAX -K_START ; k++){
		for(i = 0; i< input_E->n_conf ; i++){
		 matrix_cluster[k].a[i] = input_corr[k].a[i]*exp(Nx/2.0*(input_E[k].a[i]))/(cosh((Nx/2.0-(k+K_START))*input_E[k].a[i]))/2.0;
		 matrix_cluster[k].a[i] = sqrt(matrix_cluster[k].a[i]);
		}
	matrix_cluster[k].mean  = sqrt(exp(Nx/2.0*(input_E[k].mean))/(cosh((Nx/2.0-(k+K_START))*input_E[k].mean))/2.0);
	}


}
#endif

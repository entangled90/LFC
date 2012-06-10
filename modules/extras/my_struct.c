#include "struct.h"

#ifndef MY_STRUCT_C
#define MY_STRUCT_C

#include <stdlib.h>
#include <stdio.h>
#include <varie.h>
extern void init_int_var( rtn_int_var rtn){
	rtn.int_flat = 0.0;
	rtn.var_flat = 0.0;
	rtn.int_root = 0.0;
	rtn.var_root = 0.0;
	rtn.var_gauss = 0.0;
	rtn.int_gauss = 0.0;
	}
/* Inizializza un array di cluster, allocando la memoria necessaria */
extern void init_cluster_jk ( cluster_jk * c, int n_cluster , int vector_length){
	int i;
	for( i = 0 ;  i < n_cluster ; i++){
		c[i].n_conf = vector_length;
		c[i].a = malloc(vector_length*sizeof(double));
		c[i].error = 0.0;
	}
	}

extern void free_cluster_jk ( cluster_jk *c , int length){
	int i;
	for( i = 0; i< length; i++){
		free( (c+i)->a);
		}
	free(c);
	}

extern double variance_cluster_jk ( cluster_jk *c){
	int i ;
	double tmp = 0;
	for(i= 0 ; i < c->n_conf; i++){
		tmp +=(c->a[i] - c->mean)*(c->a[i] - c->mean);
	}
	tmp *=( c->n_conf-1);
	tmp /= c->n_conf;
	return (tmp);
	}

	/* NB: la lunghezza di vector deve essere uguale a input->n_conf */
extern void clusterize( cluster_jk *input , double *vector){
		int i = 0;
		input->mean =meanOfDoubleArray( vector, input->n_conf);
		for ( i = 0 ; i < input->n_conf ; i++){
			(input->a)[i] = input->mean + ( vector[i]-input->mean )/(double)(input->n_conf-1);
		}
}
	


#endif

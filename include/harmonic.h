#ifndef HARMONIC_H
#define HARMONIC_H

#ifndef HARMONIC_C

#include "struct.h"
extern double elagrangian(double,double);
extern double potential (double x );extern double efirst_action (double *x, int n);
extern double edelta_action ( double *x, double x_new,int position);

extern void metropolis( double *x , int position, double *x_new);
extern double correlation (double x[],int dK);
extern void DeltaE_cluster (  cluster_jk * correlation, cluster_jk *energy);
extern void matrix_element_cluster ( cluster_jk * input ,cluster_jk * input_cor, cluster_jk *matrix_cluster);
#endif
#endif

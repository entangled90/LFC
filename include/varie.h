#ifndef VARIE_H
#define VARIE_H

#ifndef VARIE_C
#include "struct.h"
extern float meanOfFloatArray(float*, int);
extern double meanOfDoubleArray (double*,int);
extern double flatPdf(double , double, double);
extern double gaussPdf(double, double,double);
extern double gaussMomentPdf(double , double , int * , double);
extern double root_exp_pdf(double x );

extern void fprintStruct ( rtn_int_var input[],int nStruct, const char* name );
extern void fprintPlot ( rtn_int_var input[],int, int ,const char* name);

extern noise fitNoise (rtn_int_var input);
extern void fprintNoisePlot ( noise input[],int Ncycle, int method, const char* name);
extern void vec_init( double *, int);
extern void clusterize( cluster_jk * , double *);
extern void vector_copy (double * vector_input, double * vector_output, int length);
#endif

#endif

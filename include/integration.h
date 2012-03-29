#ifndef INTEGRATION_H
#define INTEGRATION_H
#ifndef INTEGRATION_C


extern double Simpson ( double , double , double (*) (double) );
extern double trapezio ( double, double , double (*) (double) );

extern double partition (double, double, int, int, double (*) (double ) );


extern double gaussianQuad ( double min , double max , double (*f) (double));
#endif

#ifndef MONTECARLO_C

#include "struct.h"
extern rtn_int_var campionamentoImportanza (  double , double , int , double (*) (double ));
#endif
#endif

#define PDF_C
#include "random.h"
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589793



void exp_dble(double rd[],int n){
  int k = 0;
  ranlxd(rd,n);
  for ( k=0; k<n;k++){
    rd[k] = -log(1-rd[k]);
  }
  
}

void root_exp_dble(double rd[],int n){
   int k;
   double *ud1;
   double *ud2;
   double x1,x2,y1;
   ud1 = malloc(n*sizeof(double));
   ud2 = malloc(n*sizeof(double));
   for (k=0;k<n;k++)
   {
      gauss_dble(ud1,n);
      exp_dble( ud2,n);
      x1=ud1[k];
      x2=ud2[k];;
      y1= x1*x1+x2;
      rd[k] = y1;
     }
    free(ud1);
    free(ud2);
}


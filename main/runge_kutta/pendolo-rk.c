#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define H 0.001
/* Friction coefficient */
#define Q 0.3
/* Extern force amplitude */
#define B 1.4
#define OMEGA_EXT 0.6667
#define N_STEPS 100000000


void kn_fill (double *k1_vec, double *k2_vec,double x1,double x2, double t , double (*f1) (double ,double,double ), double (*f2) (double, double,double)){
  int i = 0;
  k1_vec[0]= H*f1(x1,x2,t);
  k2_vec[0]= H*f2(x1,x2,t);
  for(i = 0; i< 2 ; i++){
    k1_vec[i+1]= H*f1(x1+k1_vec[i]/2,x2,t+H/2 ); 
    k2_vec[i+1]= H*f2(x1,x2+k2_vec[i]/2,t+H/2 ); 
  }
  k1_vec[3]= H * f1(x1+k1_vec[2],x2,t+H);
  k2_vec[3]= H * f2(x1,x2+k1_vec[2],t+H);

}

double f2 ( double x1, double x2, double time){
  return ( -sin(x1) -Q*x2+B*cos(OMEGA_EXT*time));
}
double f1 ( double x1, double x2, double time){
  return (x2);

}
double increment_k ( double * k ){
  return ( k[0]/3+k[1]/6+k[2]/6+k[3]/3);
}

int main (int argc, char* argv[]){
  int i = 0;
  double t =0;
  double x1,x2;
  double k1[4],k2[4];
  x1 = 0;
  x2= 0;
  double tmp1,tmp2;
FILE *fp = fopen("../data/differential_equation/runge_kutta.dat","w");
for( i = 0; i< N_STEPS ; i++){
    kn_fill(k1,k2,x1,x2,t+i*H,f1,f2);
     
    tmp1 =increment_k(k1);
    tmp2 = increment_k(k2);
    if( i%10 == 0)
    fprintf(fp,"%e \t %e\n", x1,x2);
       x1 += tmp1;
       x2 += tmp2;
    
   } 
      fclose(fp);
      return(EXIT_SUCCESS);

}

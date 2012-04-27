#include <stdio.h>
#include <stdlib.h>
#include <types.h>
#include <gsl/gsl_matrix_complex_double.h>


#define R_MIN 0.4
#define R_MAX 0.5
#define V_MAX 0.2


double V_step_tunnel (double x , double y ){
	if ( ( x*x + y*y > R_MIN*R_MIN ) && (x*x +y*y < R_MAX*R_MAX))
		return V_MAX;
	else
		return 0 ;
	}
	
int main (int argc, char *argv[]){
	gsl_matrix_complex *psi;
	psi = gsl_matrix_complex_alloc (500,500);
	gsl_matrix_complex_free(psi);
	return(EXIT_SUCCESS);
	}

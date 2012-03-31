#define MONTECARLO_C
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"
#include "varie.h"
#include "struct.h"
#include "integration.h"
/*
 *
 * N è il numero di punti sparati con ogni montecarlo.
 * DEVE essere un multiplo di nRandom , per evitare errori dati dalla divisione
 * N/nRandom
 */


rtn_int_var campionamentoImportanza (  double a, double b , int N ,double (*f) (double)){
	int nRandom=50;
	int j=0;
	double *rd;
	rd = malloc(nRandom*sizeof(double));
	int i = 0;
	double integral_true = partition(a,b,N,1,f);
	rtn_int_var rtn ;
	double tmp = 0.0;
	/*
	 *	1 -> Flat
	 *	2 -> Gauss
	 *  3 -> sqrt(x) e^(-x)
	 */
	/*
	 * Mappo i random generati in [0,1] -> [a,b]
	*/
	for(j=0;j< N/nRandom;j++){
		ranlxd(rd,nRandom);
		for (i = 0; i< nRandom ; i++){
			rd[i] = a + (b-a)*rd[i];
			tmp = ( f( rd[i] ) / flatPdf( a,b ,1 ));
			rtn.var_flat += (integral_true - tmp)*(integral_true-tmp)/(N-1)/N;
			rtn.int_flat += tmp/N;
		}
	}

	 /*
	 * GAUSSSSSS
	 */
	for(j=0;j<N/nRandom;j++){
	gauss_dble(rd,nRandom);
		for ( i = 0; i< nRandom; i++){
			tmp = ( f( rd[i] ) / gaussPdf( 0.5,0,rd[i]));
			rtn.int_gauss += tmp/(double) N;
			rtn.var_gauss += (integral_true - tmp)*(integral_true-tmp)/(double )N/(N-1);
			
		}
	}
	/*
	 *
	 *ROOOOOT
	 */
	for(j=0;j<N/nRandom;j++){
	root_exp_dble(rd,nRandom);
			for ( i = 0; i< nRandom ; i++){
				tmp= ( f( rd[i] ) / root_exp_pdf(rd[i]));
				rtn.var_root +=(integral_true - tmp)*(integral_true-tmp)/(double )N/(N-1);
				rtn.int_root += tmp /(double) N;
			}
		}
		rtn.int_root*=2;
	rtn.Npnt = N;
	free(rd);
	return (rtn);
}

#define MONTECARLO_C
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"
#include "varie.h"
#include "struct.h"
#include "integration.h"
#define Nbintrap 100000
/*
 *
 * N è il numero di punti sparati con ogni montecarlo.
 * DEVE essere un multiplo di nRandom , per evitare errori dati dalla divisione
 * N/nRandom
 */


rtn_int_var campionamentoImportanza (  double a, double b , int N ,double (*f) (double)){
	int nRandom=1000;
	int j=0;
	double *rd;
	int i = 0;
	double integral_true = partition(a,b,Nbintrap,1,f);
	rtn_int_var rtn ;
	double tmp = 0.0;
	rd = malloc(nRandom*sizeof(double));
	init_int_var(rtn);
	rtn.Npnt = N;
	/*
	 *	1 -> Flat
	 *	2 -> Gauss
	 *  3 -> sqrt(x) e^(-x)
	 */
	
	/* Il doppio ciclo serve a non usare troppa memoria. In questo modo, infatti,
	 * il vettore di numeri random ha lunghezza 100, qualsiasi sia il numero
	 * dei punti utilizzati dal montecarlo.
	 * ranlxd viene chiamata con un numero maggiore di 32, in quanto consigliato dalle librerie.
	 *  Ho scelto il valore 50 per nRandom per non avere problemi con la divisione intera
	 * visto che i punti utlilizzati dal monte carlo sono multipli di 50.
	 */
	for(j=0;j< N/nRandom;j++){
		ranlxd(rd,nRandom);
		for (i = 0; i< nRandom ; i++){
			rd[i] = a + (b-a)*rd[i];
			tmp = ( f( rd[i] ) / flatPdf( a,b ,1 ));
			rtn.var_flat += (integral_true - tmp)*(integral_true-tmp)/((double)(N-1)*(double)N);
			rtn.int_flat += tmp/(double)N;
		}
	}
	 /*
	 * GAUSSSSSS
	 */
	for(j=0;j<N/nRandom;j++){
	gauss_dble(rd,nRandom);
		for ( i = 0; i< nRandom; i++){
			tmp = ( f( rd[i] ) / gaussPdf( 0.5,0,rd[i]));
			rtn.int_gauss += tmp ;
			rtn.var_gauss += (integral_true - tmp)*(integral_true-tmp)/((double )N*(double)(N-1));
			
		}
	}
	rtn.int_gauss /= (double) N ;
	/*
	 *
	 *ROOT
	 */
	for(j=0;j<N/nRandom;j++){
	root_exp_dble(rd,nRandom);
			for ( i = 0; i< nRandom ; i++){
				tmp=  f( rd[i] ) / root_exp_pdf(rd[i]);
				rtn.var_root +=(integral_true - tmp)*(integral_true-tmp)/((double )N* (double)(N-1));
				rtn.int_root += tmp /(double) N;
			}
		}
	/* Viene raddoppiato perchè la pdf estrae solo tra 0 e +inf.
	 * Utilizzabile solo per funzioni pari
	 */
	rtn.int_root*=2;

	free(rd);
	return (rtn);
}

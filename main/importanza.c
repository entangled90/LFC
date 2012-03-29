#define MAIN_PROGRAM 

#include <stdio.h>
#include <stdlib.h>
#include "integration.h"
#include "varie.h"
#include "random.h"
#include "time.h"
#include "struct.h"
#include "math.h"
#define Nstart 10000
//#include "global.h"
/*
 * Variabili Globali
 */
	static double gaussianVariance = 0.5;
	static double gaussianMean = 0;
	static int *nMoment = 0;




	/*
	 * Funzioni
	 */
double integrand (double x ){
	return(gaussMomentPdf(gaussianVariance,gaussianMean,nMoment,x));
}

int main (){
	double min =0 ;
	double max= 0;
	double *vec;
	int N = Nstart;
	vec = malloc(N*sizeof(double));
	rtn_int_var result;
	/*FILE *fp;
	fp = fopen("../data/dati.dat","w");
	
	printf("Inserisci minimo : \n");
	scanf("%lf",&min);
	printf("inserisci massimo : \n");
	scanf("%lf",&max);
	*/
	rlxd_init(1,time(NULL));
	min=-100;
	max=100;
	printf("\t \t  \t flat  \t gauss \t trapezi \t root_exp \n");
	int i = 0;
	nMoment=&i;
	
	for ( N= Nstart; N<Nstart*10; N*=3){
	  printf("I punti usati sono %d \n",N);
	 for(i=2; i < 7; i+=2){
		result = campionamentoImportanza(min,max,N,integrand);
		printf("Il momento  numero %d : %lf \t %lf \t %lf \t \n",i,
		result.int_flat, result.int_gauss,result.int_root);
		printf("\n Errori: \t");
		printf("\t %lf \t %lf \t %lf \n",sqrt(result.var_flat),sqrt(result.var_gauss),sqrt(result.var_root));
		printf("\n");
		}
	}
	/*
	for ( i= 0; i< N ;i++){
	fprintf(fp,"%lf\n",vec[i]);
	}
	fclose(fp);
	*/
	
	free(vec);
	return(EXIT_SUCCESS);
}

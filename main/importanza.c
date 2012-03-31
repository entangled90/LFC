#define MAIN_PROGRAM 

#include <stdio.h>
#include <stdlib.h>
#include "integration.h"
#include "varie.h"
#include "random.h"
#include "time.h"
#include "struct.h"
#include "math.h"
#include "plot.h"
#define Nstart 1000
#define  Ncycle 6
#define MaxMoment 4

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

	int N = Nstart;
	rtn_int_var result[Ncycle*MaxMoment/2];
	noise noise_result[Ncycle*MaxMoment/2];
	int momentIndex = 1;

	/*
	printf("Inserisci minimo : \n");
	scanf("%lf",&min);
	printf("inserisci massimo : \n");
	scanf("%lf",&max);
	*/
	rlxd_init(1,time(NULL));
	min=-100;
	max=100;
	printf("\t \t  \t flat  \t gauss  \t root_exp \n");
	int i = 0;
	int j = 0;
	nMoment=&momentIndex;
	for ( j=0; j<Ncycle; j++){
		printf("I punti usati sono %d \n",N);
		momentIndex = 2;
		for(i=0; i < MaxMoment/2; i++){
			result[(MaxMoment/2)*j+i] = campionamentoImportanza(min,max,N,integrand);
			result[(MaxMoment/2)*j+i].Npnt = N;
			printf("Il momento  numero %d : %lf \t %lf \t %lf \t \n",momentIndex,
			result[(MaxMoment/2)*j+i].int_flat, result[(MaxMoment/2)+i].int_gauss,result[(MaxMoment/2)*j+i].int_root);
			printf("\n Errori: \t");
			printf("\t %lf \t %lf \t %lf \n",sqrt(result[(MaxMoment/2)*j+i].var_flat),sqrt(result[(MaxMoment/2)*j+i].var_gauss),sqrt(result[(MaxMoment/2)*j+i].var_root));
			printf("\n");
			momentIndex+=2;
			}
		N*=2;
	}
	for(j=0; j< Ncycle*MaxMoment/2; j++){
	noise_result[j] = fitNoise(result[j]);
	}
	/*
	 *Stampano vari file:
	 * plot contiene tutti i risultati
	 * printplot per il valore del'integrale con errore
	 * fprintnoiseplot per il valore del rumore
	 */
	fprintStruct(result,Ncycle*MaxMoment/2,"../data/plot.dat");
	fprintPlot(result,Ncycle*MaxMoment/2,1,"../data/flat.dat");
	fprintPlot(result,Ncycle*MaxMoment/2,2,"../data/gauss.dat");
	fprintPlot(result,Ncycle*MaxMoment/2,3,"../data/root.dat");
	fprintNoisePlot(noise_result,Ncycle*MaxMoment/2,1,"../data/flat_scaled.dat");
	fprintNoisePlot(noise_result,Ncycle*MaxMoment/2,2,"../data/gauss_scaled.dat");
	fprintNoisePlot(noise_result,Ncycle*MaxMoment/2,3,"../data/root_scaled.dat");
	plot("../data/flat.dat","../data/gauss.dat","../data/root.dat");
	plot_noise("../data/flat_scaled.dat","../data/gauss_scaled.dat","../data/root_scaled.dat");
	return(EXIT_SUCCESS);
}

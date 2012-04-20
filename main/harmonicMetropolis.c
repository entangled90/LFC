#define MAIN_PROGRAM
#include <stdio.h>
#include <stdlib.h>
#include "harmonic.h"
#include <sys/ioctl.h>
#include "random.h"
#include "time.h"
#include "plot.h"
#include "globalharmonic.h"
#include "varie.h"
#include "math.h"
#include "struct.h"

int main(){
	int i = 0;
	int j = 0;
	int k = 0;
	double x[Nx];
	double xnew[Nx];
	double O[Nx];
	double *mean_bin;
	double *variance_bin;
	double *deltaE_vector;
	double temp;
	struct winsize w;
	cluster_jk *correlation_clustered;
	cluster_jk *energy_clustered;
	/*
	 * Check per non avere problemi con la divisione intera
	*/
	if(N_BIN_CHECK != 0){
		printf("N_BIN_CHECK ERROR!");
		return(EXIT_FAILURE);
	}
	/* Necessario per importare il numero di colonne della finestra di terminale*/
	ioctl(0,TIOCGWINSZ,&w);
	/* Allocazione della memoria*/
	mean_bin = malloc(N_BIN*(Nx)*sizeof(double));
	variance_bin = malloc((Nx)*sizeof(double));
	deltaE_vector = malloc((K_MAX-K_START)*sizeof(double));
	correlation_clustered = malloc(K_MAX*sizeof(cluster_jk));
	energy_clustered = malloc((K_MAX-1)*sizeof(cluster_jk));
	/* Alloca i vettori per la struct cluster_jk e impone che il numero di configurazioni è N_BIN*/
	init_cluster_jk(correlation_clustered , K_MAX , N_BIN);
	init_cluster_jk(energy_clustered, K_MAX -2, N_BIN);
	/*
	 * Inizializzazione dei vettori a zero
	 */
	vec_init(x,Nx);
	vec_init(O,Nx);
	vec_init(mean_bin,N_BIN*(Nx));
	vec_init(variance_bin, Nx);
	/*
	 * Init random generator
	 */
	rlxd_init(1,time(NULL));
	printf("[");
	/* Ciclo sugli sweep */
	for ( j = 0; j< N_SWEEP; j++){
		ranlxd(xnew,Nx);
		/* Esegue uno sweep completo su tutte le coordinate */
		for ( i = 0 ; i< Nx; i++){
			xnew[i] = x[i]+ 2*DELTA*(xnew[i]-0.5);
			metropolis(x,i,xnew[i]);
		}
		/* Calcola la media delle correlazione su ogni sweep dopo i 200 */
		if ( j > THERM_CONST){
			/*Calcola i valori delle correlazioni, la media per ogni bin, su ogni valore di |l-k|*/
			for(k = 0; k< Nx ; k++){
				temp =  correlation (x,k);
				O[k] += temp/(double) (N_SWEEP-THERM_CONST) ;
				mean_bin[((j-THERM_CONST)/BIN_WIDTH) + k*N_BIN] += temp/(double) BIN_WIDTH ;
			}
		}
	/* Stampa a terminale la progressione del programma tramite una barra  */
			if ( (int) ((j*1.0)/N_SWEEP*(w.ws_col-1))  != (int) ((j-1)*1.0/N_SWEEP*(w.ws_col-1)) ){
				printf("=");
				fflush(stdout);
			}
	}
	printf("]\n");
	/* Calcola la deviazione standard della media */
	for ( i = 0 ; i < Nx; i++){
	
		for( j  = 0 ; j< N_BIN ; j++){
			variance_bin[i] += pow((mean_bin[j+i*N_BIN] -O[i]),2.0) / (double) (N_BIN-1)/(double) N_BIN; 
		}
		variance_bin[i] = sqrt(variance_bin[i]);
	}
	for ( i = 0 ; i < K_MAX ; i++){
		vector_copy( &mean_bin[(i+1)*N_BIN], correlation_clustered[i].a,N_BIN);
		clusterize( (correlation_clustered+i) );
	}
	 DeltaE_cluster ( correlation_clustered , energy_clustered, deltaE_vector);
	 //printf("%lf \n " , correlation_clustered[2].mean);
	{
	/*
	 *	Calcola il valore   dell'energia Delta E = E_1 - E_0 per le prime quattro correlazioni,
	 * le uniche per cui si ha segnale.
	 *	Salva i valori nel file contenente.
	 */
	 for(i = 0; i< K_MAX - K_START; i++)
		clusterize((energy_clustered+i));
	
	FILE *fp_energy = fopen("../data/harmonic/energy.dat","w");
	for ( i = 0 ; i< K_MAX-K_START ; i++){
		fprintf(fp_energy," %d \t %14.10e \n", i+1 , deltaE_vector[i]);
	}
	/* File che contiene i valori della correlazione per i valori | l -k| */
	FILE *fp = fopen("../data/harmonic/harmonic.dat","w");
	for(i = 0; i<Nx; i++){
		fprintf(fp,"%d \t %lf \n", i, O[i]);
	}
	/* FIle con la varianza per ogni osservabile di correlazione */
	FILE *fp_var = fopen("../data/harmonic/harm_var_bin.dat","w");
	for ( i = 0 ; i < Nx; i++)
		fprintf(fp_var,"%d \t %14.10e \n", i , variance_bin[i]);
	fclose(fp);
	fclose(fp_var);
	plot_harmonic("../data/harmonic/harmonic.dat","../data/harmonic/harmonic.eps");
	return(EXIT_SUCCESS);
	}
}

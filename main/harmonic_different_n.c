#define MAIN_PROGRAM
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include "harmonic.h"
#include "random.h"
#include "time.h"
#include "plot.h"
#include "globalharmonic.h"
#include "varie.h"
#include "math.h"
#include "struct.h"

int main(int argc, char *argv[]){
	int i = 0;
	int j = 0;
	int k = 0;
	int N_SWEEP = 10000;
	double x[Nx];
	double xnew[Nx];
	double O[Nx];
	double *mean_bin;
	double *variance_bin;
	double temp;
	struct winsize w;
	double deltaE_variance = 0.0;
	double matrix_element_variance = 0.0;
	int N_BIN;

	cluster_jk *correlation_clustered;
	cluster_jk *energy_clustered;
	cluster_jk *energy_time_clustered;
	cluster_jk *matrix_clustered;
	cluster_jk *matrix_time_clustered;
	/*
	 * Check per non avere problemi con la divisione intera
	*/
	if( argc == 1){
		printf("Parameter Error\n");
		return(EXIT_FAILURE);
	}
	else {
		N_SWEEP = atoi(argv[1]); 
		N_BIN =  ((N_SWEEP)/BIN_WIDTH);
	}
	if( N_SWEEP % BIN_WIDTH != 0){
		printf("N_BIN_CHECK ERROR!\n");
		return(EXIT_FAILURE);
	}
	/* Necessario per importare il numero di colonne della finestra di terminale*/
		ioctl(0,TIOCGWINSZ,&w);
	/* Allocazione della memoria*/
		mean_bin = malloc(N_BIN*(Nx)*sizeof(double));
		variance_bin = malloc((Nx)*sizeof(double));
		//deltaE_vector = malloc((K_MAX-K_START)*sizeof(double));
		correlation_clustered = malloc((K_MAX-K_START + 2)*sizeof(cluster_jk));
		energy_clustered = malloc((K_MAX-K_START)*sizeof(cluster_jk));
		energy_time_clustered = malloc(sizeof(cluster_jk));
		matrix_clustered = malloc((K_MAX-K_START)*sizeof(cluster_jk));
		matrix_time_clustered = malloc(sizeof(cluster_jk));
		
	/* Alloca i vettori per la struct cluster_jk e impone che il numero di configurazioni sia N_BIN*/
		init_cluster_jk(correlation_clustered , K_MAX , N_BIN);
		init_cluster_jk(energy_clustered, K_MAX -K_START, N_BIN);
		init_cluster_jk(energy_time_clustered, 1 , N_BIN);
		init_cluster_jk(matrix_clustered,K_MAX - K_START , N_BIN);
		init_cluster_jk(matrix_time_clustered, 1 , N_BIN);
	/* Inizializzazione dei vettori a zero */
		vec_init(x,Nx);
		vec_init(O,Nx);
		vec_init(mean_bin,N_BIN*(Nx));
		vec_init(variance_bin, Nx);
	/*Init random generator */
		rlxd_init(1,time(NULL));
		printf("Metropolis' Progression: \n[");
	/* Prepara il vettore, facendo THERM_CONST sweep prima di iniziare a calcolare */{
		
		for ( j = 0; j< THERM_CONST; j++){
			ranlxd(xnew,Nx);
		/* Esegue uno sweep completo su tutte le coordinate */
			for ( i = 0 ; i< Nx; i++){
				xnew[i] = x[i]+ 2*DELTA*(xnew[i]-0.5);
				metropolis(x,i,xnew+i);
			}
		}
		for ( j = 0; j< N_SWEEP; j++){
			ranlxd(xnew,Nx);
		/* Esegue uno sweep completo su tutte le coordinate */
			for ( i = 0 ; i< Nx; i++){
				xnew[i] = x[i]+ 2*DELTA*(xnew[i]-0.5);
				metropolis(x,i,xnew+i);
			}
			/*Calcola i valori delle correlazioni, la media per ogni bin, su ogni valore di |l-k|*/
				for(k = 0; k< Nx ; k++){
					temp =  correlation (x,k);
					O[k] += temp/(double) (N_SWEEP-THERM_CONST) ;
					mean_bin[((j)/BIN_WIDTH) + k*N_BIN] += temp/(double) BIN_WIDTH ;
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
	/* Copia il vettore delle correlazioni binnato nei cluster e clusterizza le strutture */
		for ( i = 0 ; i < K_MAX ; i++){
			clusterize( (correlation_clustered+i),&mean_bin[(i+1)*N_BIN]  ); // clusterize calcola anche la media
		}
	/*	Calcola DeltaE per ogni configurazione e in valor medio.
	 * 	Il valor medio (per ogni tempo) viene salvato nella media della struttura
	 */
		DeltaE_cluster ( correlation_clustered , energy_clustered);
		matrix_element_cluster(energy_clustered, correlation_clustered, matrix_clustered);
	/*
	 *	Calcola il valore   dell'energia Delta E = E_1 - E_0 per le prime quattro correlazioni,
	 * le uniche per cui si ha segnale.
	 *	Salva i valori nel file contenente.
	 */
		/*Inizializza il vettore della struct energy_time_clustered a 0 */
		vec_init(energy_time_clustered->a, energy_time_clustered->n_conf);
		vec_init(matrix_time_clustered->a, matrix_time_clustered->n_conf);
		/* energy_time_clustered contiene la media tra i tempi markoviani di DeltaE^j*/
		energy_time_clustered->mean = 0.0;
		matrix_time_clustered->mean = 0.0;
		for(i = 0; i< K_MAX - K_START; i++){
			for( j = 0; j< N_BIN; j++){
				energy_time_clustered->a[j] += energy_clustered[i].a[j]/(double)(K_MAX-K_START);
				matrix_time_clustered->a[j] += matrix_clustered[i].a[j]/(double)(K_MAX-K_START);
			}
			energy_time_clustered->mean += energy_clustered[i].mean/(double) (K_MAX -K_START);
			matrix_time_clustered->mean += matrix_clustered[i].mean/(double) (K_MAX -K_START);
			printf("Elemento clustered : %e \n", matrix_clustered[i].mean);
		}
		deltaE_variance = variance_cluster_jk(energy_time_clustered );
		matrix_element_variance = variance_cluster_jk ( matrix_time_clustered);
		
	/* Stampa su file i valori medi dell'energia 
	*/
		FILE *fp_energy = fopen("../data/harmonic/n_variance_energy.dat","a");
		FILE *fp_matrix = fopen("../data/harmonic/n_variance_matrix.dat","a");
			fprintf(fp_energy,"%d \t %14.10e \n",N_SWEEP, sqrt(deltaE_variance));
			fprintf(fp_matrix,"%d \t %14.10e \n", N_SWEEP, sqrt(matrix_element_variance));
	/*	for ( i = 0 ; i< K_MAX-K_START ; i++){
			printf("Energy clustered : %14.10e \n", energy_clustered[i].mean);
		}
	*/
		printf("Delta E medio: \t %14.10e \t Dev. std: \t %e \n", energy_time_clustered->mean, sqrt(deltaE_variance));
		printf("Elemento di matrice medio: %14.10e\tDev. std:%e \n", matrix_time_clustered->mean, sqrt(matrix_element_variance));
	/*
		FILE *fp_matrix_element = fopen("../data/harmonic/matrix_element.dat","a");
		FILE *fp_matrix_variance = fopen("../data/harmonic/matrix_variance.dat","a");
		fprintf(fp_matrix_element,"%14.10e\n", matrix_time_clustered->mean);
		fprintf (fp_matrix_variance,"%14.10e\n", sqrt(matrix_element_variance));
		
	/* File che contiene i valori della correlazione per i valori | l -k| 
		FILE *fp = fopen("../data/harmonic/harmonic.dat","w");
		for(i = 0; i<Nx; i++){
			fprintf(fp,"%d \t %lf \n", i, O[i]);
		}
	/* FIle con la varianza per ogni osservabile di correlazione 
		FILE *fp_var = fopen("../data/harmonic/harm_dev_std_bin.dat","w");
		for ( i = 0 ; i < Nx; i++)
			fprintf(fp_var,"%d \t %14.10e \n", i , variance_bin[i]);
	*/
		free(mean_bin);
		free(variance_bin);
	/* Funzione che fa il free del vettore contenuto nella struct e poi della struttura */
		free_cluster_jk(correlation_clustered,K_MAX);
		free_cluster_jk(energy_clustered, K_MAX - K_START);
		free_cluster_jk(energy_time_clustered ,1);
		free_cluster_jk(matrix_clustered,K_MAX - K_START);
	/*Chiude i file ed esegue il plot */
		fclose(fp_matrix);
		//fclose(fp_var);
		fclose(fp_energy);
		// plot_harmonic("../data/harmonic/harmonic.dat","../data/harmonic/harmonic.eps");
		// plot_harmonic("../data/harmonic/harm_dev_std_bin.dat","../data/harmonic/harm_dev_std_bin.eps");
		// plot_harmonic("../data/harmonic/energy.dat","../data/harmonic/energy.eps");
		// fit("../data/harmonic/energy.dat", "../data/harmonic/energy_histogram.eps" ,energy_time_clustered->mean, deltaE_variance);
		// fit("../data/harmonic/energy_variance.dat", "../data/harmonic/energy_variance_histogram.eps", deltaE_variance, 1e-5);
		return(EXIT_SUCCESS);
	}

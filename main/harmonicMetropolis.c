#define MAIN_PROGRAM
#include <stdio.h>
#include "harmonic.h"
#include <stdlib.h>
#include "random.h"
#include "time.h"
#include "plot.h"
#include "globalharmonic.h"
#include "varie.h"
#include "math.h"
#include <sys/ioctl.h>

/* Nota bene che N_BIN DEVE essere intero.
 * ħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħ
 * ħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħħ
 */
#define N_BIN (N_SWEEP-THERM_CONST)/BIN_WIDTH
#define K_MAX 4



int main(){
	int i = 0;
	int j = 0;
	int k = 0;
	//int l = 0;
	double x[Nx];
	double xnew[Nx];
	//double *action_cold;
	double O[Nx];
	double O2[Nx];
	double *mean_bin;
	double *variance_bin;
	double *deltaE_vector;
	double temp;
	struct winsize w;
	
	//double *mean;
	// action_cold = malloc(N_SWEEP*sizeof(double));
	ioctl(0,TIOCGWINSZ,&w);
	mean_bin = malloc(N_BIN*(Nx)*sizeof(double));
	//mean = malloc( (Nx)*sizeof(double));
	variance_bin = malloc((Nx)*sizeof(double));
	deltaE_vector = malloc((K_MAX-1)*sizeof(double));
	srand(time(NULL));
	
	rlxd_init(1,time(NULL));
	/* Vector initialisation to zero */
	vec_init(x,Nx);
	vec_init(O,Nx);
	vec_init(mean_bin,N_BIN*(Nx));
	vec_init(variance_bin, Nx);
	printf("[");
	for ( j = 0; j< N_SWEEP; j++){
		ranlxd(xnew,Nx);
		/*Sweep completo */
		for ( i = 0 ; i< Nx; i++){
			xnew[i] = x[i]+ 2*DELTA*(xnew[i]-0.5);
			metropolis(x,i,xnew[i]);
		}if ( j > 200){
			/* Calcola la media delle correlazione su ogni sweep dopo i 200 */
			for(k = 0; k< Nx ; k++){
				temp =  correlation (x,k)/(double) (N_SWEEP-THERM_CONST);
				O[k] += temp ;
				O2[k] += temp*temp*(double) (N_SWEEP-THERM_CONST);
				mean_bin[((j-THERM_CONST)%BIN_WIDTH)*Nx + k] += temp/(double) N_BIN;
			}	
		}
	//action_cold[j] = efirst_action(x,Nx);
	if ( (int) ((j*1.0)/N_SWEEP*(w.ws_col-1))  != (int) ((j-1)*1.0/N_SWEEP*(w.ws_col-1)) )
		printf("=");
		fflush(stdout);
	}
	printf("]\n");
	for ( i = 0 ; i < Nx; i++){
		for( j  = 0 ; j< N_BIN ; j++){
			variance_bin[i] += pow((mean_bin[j*Nx+i]*mean_bin[j*Nx+i] -O[i]),2) / (double) (N_BIN-1)/(double) N_BIN; 
		}
			variance_bin[i] = sqrt(variance_bin[i]);
	}
	DeltaE(O, K_MAX, deltaE_vector);
	/* File contenente i valori dell'energia Delta E = E_1 - E_0 per le prime quattro correlazione,
	 * le uniche per cui si ha segnale*/
	FILE *fp_energy = fopen("../data/harmonic/energy.dat","w");
	for ( i = 0 ; i< K_MAX-1 ; i++){
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
	//free(action_cold);
	plot_harmonic("../data/harmonic/harmonic.dat","../data/harmonic/harmonic.eps");
	return(EXIT_SUCCESS);
	}

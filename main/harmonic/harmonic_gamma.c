#define MAIN_PROGRAM
#include <stdio.h>
#include "harmonic.h"
#include <stdlib.h>
#include "random.h"
#include "time.h"
#include "plot.h"
#include "globalharmonic.h"

#define CORR_SIZE 50
#define N_SWEEP 10000
int main(){
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int length = N_SWEEP*Nx;
	double *x;
	double *xnew;
	double O[Nx-1];
	double O2[Nx-1];
	double O_mix[Nx-1][CORR_SIZE];
	double GammaT[Nx-1][CORR_SIZE];
	rlxd_init(1,time(NULL)%100000);
	FILE *fp = fopen("../../data/harmonic/action.dat","w");
	printf("Inserire \"1\" per azione fredda, \"2\" per azione calda\n");
	scanf("%d",&l);
	/*initialisation */
	x=  malloc(length*sizeof(double));
	xnew =  malloc(Nx*sizeof(double));
	/*Inizializzazione a zero degli array */
	for ( i = 0; i< Nx-1; i++){
		O[i] = 0.0;
		O2[i] = 0.0;
		for(j = 0; j<CORR_SIZE ; j++){
			O_mix[i][j] = 0.0;
		}
	}
	switch(l){
	  case 1:
	    for( i = 0; i< Nx; i++)
	      x[i] = 0.0;
	    break;
	  case 2:
	    ranlxd(x,Nx);
	    for( i = 0; i< Nx; i++)
	      x[i] = (2*(x[i]-0.5))*500;
	    break;
	  default:
	    printf("Valore inserito non corretto\n");
	    exit(EXIT_FAILURE);
	}
	for ( j = 0; j< N_SWEEP; j++){
		ranlxd(xnew,Nx);
		/*Sweep completo */
		for ( i = 0 ; i< Nx; i++){
			xnew[i] = x[j*Nx+i]+ 2*DELTA*(xnew[i]-0.5);
			metropolis(x +j*Nx,i,xnew+i);
		}
		if( j < 2000)
			  fprintf(fp,"%d \t %e\n", j , efirst_action(x+j*Nx,Nx));

		if( j != N_SWEEP-1){
			for (i = 0 ; i<Nx;i++){
				*(x+(j+1)*Nx+i) =*(x+j*Nx+i);
			}
		}
		if ( j> THERM_CONST){
			/* Calcola la media delle correlazione su ogni sweep dopo i 200 */
			for(k = 1; k< Nx ; k++){
				O[k-1] += correlation (&x[j*Nx],k )/(double) (N_SWEEP-THERM_CONST);
				O2[k-1] += correlation (&x[j*Nx],k)*correlation (&x[j*Nx],k)/(double) (N_SWEEP-THERM_CONST);
				/*Calcola < O_i * O_(i-l)> per ogni valore di k*/
				for( l = 0 ; l< CORR_SIZE ; l++){
					O_mix[k-1][l] += correlation (&x[j*Nx],k)*correlation (&x[(j-l)*Nx],k )/(double) (N_SWEEP-THERM_CONST);
				}
			}
		}
		
	}
	fclose(fp);
	for ( l = 0; l<CORR_SIZE ; l++){
			for (k = 0; k<Nx-1;k++){
			GammaT[k][l] = (O_mix[k][l] - O[k]*O[k])/(-O[k]*O[k]+O2[k]);
			}
		}
	fp = fopen("../../data/harmonic/gamma_t.dat","w");
		for ( j = 0 ; j< Nx-1; j++){
			for ( i = 0 ; i < CORR_SIZE-1 ; i++){
				fprintf(fp," %d \t %lf \n", i+1, GammaT[j][i]);
			}
		}
	FILE *fp1 = fopen("../../data/harmonic/harmonic2.dat","w");
	for(i = 0; i<Nx-1; i++){
		fprintf(fp1,"%d \t %lf \n", i+1, O[i]);
	}
	fclose(fp);
	fclose(fp1);
	free(x);
	free(xnew);
	plot_harmonic("../../data/harmonic/harmonic2.dat","../../data/harmonic/harmonic2.eps");
	plot_harmonic("../../data/harmonic/gamma_t.dat","../../data/harmonic/gamma_t.eps");
	return(EXIT_SUCCESS);
	}

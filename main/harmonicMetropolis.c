#define MAIN_PROGRAM
#include <stdio.h>
#include "harmonic.h"
#include <stdlib.h>
#include "random.h"
#include "time.h"
#include "plot.h"
#include "globalharmonic.h"

#define Nsweep 2000

int main(){
	int i = 0;
	int j = 0;
	int k = 0;
	double x_cold[Nx];
	double xnew[Nx];
	double action_cold[Nsweep];
	double dKcorr[Nx];
	srand(time(NULL));
	rlxd_init(1,time(NULL));
	/* Cold initialisation */
	for ( i = 0; i< Nx; i++){
		x_cold[i]= 0.0;
		dKcorr[i] = 0.0;
	}
	for ( j = 0; j< Nsweep; j++){
		ranlxd(xnew,Nx);
		for ( i = 0 ; i< Nx; i++){
			xnew[i] = x_cold[i]+ 2*DELTA*(xnew[i]-0.5);
			metropolis(x_cold,i,xnew[i]);
		}if ( Nsweep > 200){
			for(k = 0; k< Nx ; k++){
			dKcorr[k] += correlation (x,k );
			}
		}
	action_cold[j] = efirst_action(x_cold,Nx);
	}
	FILE *fp_cold = fopen("../data/harmonic/harmonic.dat","w");
	for(i = 0; i<Nsweep; i++){
		fprintf(fp_cold,"%d \t %lf \n", i, action_cold[i]);
	}
	//plot_harmonic("../data/harmonic/harmonictestcold.dat","../data/harmonic/harmonictesthot.dat");
	return(EXIT_SUCCESS);
	}

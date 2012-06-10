#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "integration.h"
#include "varie.h"
#include "plot.h"

#define MAIN_PROGRAM
#define N_BIN_DEFAULT (100)
#define X_MIN_DEFAULT 1.0
#define X_MAX_DEFAULT 2.0
/*Funzioni da integrare */
double integrandLog( double x ){
	return ( log(1+x));
	}
double integrandPoly (double x ){
	return (pow(x,9) - pow(x,7) + 3.0);
	}

/*Primitive Esatte delle funzioni da integrare*/
double primitiveLog (double x){
	return	(-x + log(1 + x) + x*log(1 + x)) ;	
	}
double primitivePoly (double x ){
		return	(pow(x,10)/10 - pow(x,8)/8 + 3*x) ;	
	}
	
/*Valuta il valore dell'integrale, data la primitiva*/
double realIntegration (double (*primitive) (double x),double xmin, double xmax){
	return (primitive(xmax) - primitive(xmin));	
	}

int main (int argc, char *argv[]) {
	double xMin ; 
	double xMax ;
	int n ;
	int choice = 0;
	double maxf1,maxf4;
	/* Necessari per il corretto funzionamento del programma*/
	if ( argc < 4 ){
		n = N_BIN_DEFAULT;
		printf("ATTENZIONE: il programma va lanciato con la sintassi: integral xMin xMax N_Bin\n");
		printf("Numero di bin impostato a %d.\nEstremo inferiore dell'integrale impostato a %e.\nEstremo inferiore dell'integrale impostato a %e.\n", N_BIN_DEFAULT, X_MIN_DEFAULT,X_MAX_DEFAULT);
		xMin=X_MIN_DEFAULT;
		xMax=X_MAX_DEFAULT;
	}
	else{
		xMin=  atof(argv[1]);
		xMax = atof(argv[2]);
		n = atoi(argv[3]);
	}
	double (* integrand ) (double  );
	double (* primitive) (double );
	printf("Scegli la funzione da integrare : (1) = Log , (2) = Poly \n");
	scanf("%d", &choice);
	switch(choice) {
		case 1 :
			integrand = (integrandLog);	 //, double xMin, double xMax);
			primitive = primitiveLog;
			maxf1= 0.5;
			maxf4 = -0.0740742;
			
		break;
		case 2 : 
			integrand = (integrandPoly);
			primitive = primitivePoly;
			maxf1 = 1856;
			maxf4 = 90048;
		break;		
		default:
			printf("Inserisci 1 o 2 ");
			return (EXIT_FAILURE);
			break;
	}
	/* Calcolo integrali con il numero di punti inserito in tutti e 3 i modi */
	double integral_trap = partition(xMin,xMax, n, 1, integrand );
	double integral_simp = partition(xMin,xMax, n, 2, integrand );
	double integral_quad = partition(xMin,xMax, n,3,integrand);
	double integral_true = realIntegration(primitive , xMin,xMax);
	printf("Valore esatto \t Trapezi \t Simpson \t gaussQuad\n");
	printf("%e \t %e \t %e \t %e\n", integral_true , integral_trap, integral_simp,integral_quad );
	printf("Err. Trap \t Err. Simp \t Err.Gauss\n%e \t %e \t %e \n ", (integral_true - integral_trap) , (integral_true - integral_simp), integral_true-integral_quad );
	printf("\n");
/* Nel caso viene passato un'ulteriore argomento al programma, esso ciclerà dal numero di punti inserito come 3° argomento, fino al
 * numero di punti inserito come 4° argomento.
 * Necessario per valutare l'andamento asintotico degli algoritmi
 */
	if(argc == 5){
		int n_loops = atoi(argv[4]);
		FILE *fp_trap = fopen("../../data/integral/trap.dat","w");
		FILE *fp_simp = fopen("../../data/integral/simp.dat","w");
		FILE *fp_gauss = fopen("../../data/integral/gauss.dat","w");
		FILE *fp_true = fopen("../../data/integral/true.dat","w");
		FILE *fp_total = fopen("../../data/integral/total.dat","w");
		int i = 0;
		for ( i = n; i< n_loops; i*=2 ){
			fprintf(fp_trap,"%d\t%e\t%e\n",i,fabs(integral_true - partition(xMin,xMax, i, 1, integrand)),fabs(pow((xMax-xMin)/i,2)*(xMax-xMin)/12*maxf1));
			fprintf(fp_simp,"%d\t%e\t%e\n",i,fabs(integral_true -partition(xMin,xMax, i, 2, integrand)),fabs(pow((xMax-xMin)/i,4)*(xMax-xMin)/90*maxf4));
			fprintf(fp_gauss,"%d\t%e\n",i,integral_true -partition(xMin,xMax,i,3,integrand ));
			fprintf(fp_true,"%d\t%e\n",i,realIntegration(primitive,xMin,xMax));
			fprintf(fp_total,"%d\t%e\t%e\t%e\n",i,integral_true - partition(xMin,xMax, i, 1, integrand ),integral_true - partition(xMin,xMax, i, 2, integrand ),integral_true - partition(xMin,xMax, i, 3, integrand ));
			
			printf("Numero di intervalli: %d\n",i);
		}
		fclose(fp_trap);
		fclose(fp_gauss);
		fclose(fp_true);
		fclose(fp_simp);
	}
	//plot("../../data/integral/trap.dat","../../data/integral/simp.dat","../../data/integral/gauss.dat");
	return(EXIT_SUCCESS);
}

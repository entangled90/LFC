#define VARIE_C
#include <math.h>
#include "varie.h"
#include "struct.h"
#include "stdio.h"
#define PI  4*atan(1)

	float meanOfFloatArray( float array[], int n){
	int i = 0;
	float sum = 0;
	for (i= 0; i< n; i++){
		sum += array[i];
		
		}
	return ( sum/n);
	}


	double meanOfDoubleArray( double array[] , int n){
	int i = 0;
	double sum = 0;
	for (i= 0; i< n; i++){
		sum += array[i];
		
		}
	return ( sum/n);
	}

	double flatPdf(double min, double max,double x){
		return  1/(max-min);
	}

	double gaussPdf(double variance, double mean, double x){
		return  exp(-(mean-x)*(mean-x)/(2*variance))/sqrt(2*PI*variance) ;
	}

	double gaussMomentPdf(double variance, double mean , int * n , double x){
	return pow(x,*n)/sqrt(2*PI*variance)*exp(-(mean-x)*(mean-x)/(2*variance));
	}
	
	double root_exp_pdf(double x){
	    return 2/sqrt(PI)*exp(-x)*sqrt(x);
	}


	void fprintStruct (rtn_int_var input[],int nStruct, const char *name ){
	FILE *fp;
	fp = fopen(name,"w");
	int i=0;
	fprintf(fp,"# N \t NTFlat \t NT Gauss \t NT root\tDev Flat\tDev Gauss\t Dev Root \n");	
	for(i=0; i< nStruct; i++){
	/* Fa il printf della struct scrivendo :
	 * N 	int-flat	 int-gauss		int-root	dev-flat	dev-gauss	dev-root
	 */
		fprintf(fp,"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", input[i].Npnt,input[i].int_flat,
		input[i].int_gauss, input[i].int_root, sqrt(input[i].var_flat),sqrt(input[i].var_gauss),sqrt(input[i].var_root) );
		}
	fclose(fp);
	}

	
	 void fprintPlot ( rtn_int_var input[],int Ncycle, int method, const char* name){
	FILE *fp;
	fp = fopen(name,"w");
	double print_int=0.0;
	double print_var=0.0;
	int i=0;
	for(i=0; i< Ncycle; i++){
	/* Fa il printf della struct scrivendo :
	 * N 	int(method) dev (method)
	*/
			/*
			* 1 -> Flat
			* 2 -> Gauss
			* 3 -> Root
			*/
			switch(method){
			case 1:
				print_int =input[i].int_flat;
				print_var = input[i].var_flat;
				break;
			case 2 :
				print_int =input[i].int_gauss;
				print_var = input[i].var_gauss;
				break;
			case 3:
				print_int =input[i].int_root;
				print_var = input[i].var_root;
				break;
			default:
			printf("Error: bad method number in fprintPlot \n");
			break;
			}
			fprintf(fp,"%d \t %lf \t %lf\n", input[i].Npnt, print_int,print_var );
	}
	fclose(fp);
	}

	noise fitNoise (rtn_int_var input){
		noise out;
		out.noise_flat= sqrt(input.var_flat)/input.int_flat;
		out.noise_flat_scaled = out.noise_flat*sqrt(input.Npnt);
		out.noise_gauss= sqrt(input.var_gauss)/input.int_gauss;
		out.noise_gauss_scaled = out.noise_gauss*sqrt(input.Npnt);
		out.noise_root= sqrt(input.var_root)/input.int_root;
		out.noise_root_scaled = out.noise_root*sqrt(input.Npnt);
		out.Npnt = input.Npnt;
	return out;
	}

	 void fprintNoisePlot ( noise input[],int Ncycle, int method, const char* name){
	FILE *fp;
	fp = fopen(name,"w");
	double print_noise_scaled=0.0;
	int i=0;
	for(i=0; i< Ncycle; i++){
	/* Fa il printf della struct scrivendo :
	 * N 	int(method) dev (method)
	*/
			/*
			* 1 -> Flat
			* 2 -> Gauss
			* 3 -> Root
			*/
			switch(method){
			case 1:
				print_noise_scaled = input[i].noise_flat_scaled;
				break;
			case 2 :
				print_noise_scaled = input[i].noise_gauss_scaled;
				break;
			case 3:
				print_noise_scaled = input[i].noise_root_scaled;
				break;
			default:
			printf("Error: bad method number in fprintPlot \n");
			break;
			}
			fprintf(fp,"%d \t %lf \n", input[i].Npnt,print_noise_scaled );
	}
	fclose(fp);
	}


#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#ifndef MY_STRUCT_C
typedef struct {
	double int_flat;
	double int_gauss;
	double int_root;
	double var_flat;
	double var_gauss;
	double var_root;
	int Npnt ;
	} rtn_int_var ;

typedef struct {
	double noise_flat;
	double noise_gauss;
	double noise_root;
	double noise_flat_scaled;
	double noise_gauss_scaled;
	double noise_root_scaled;
	int Npnt ;
	} noise ;
	
	extern void init_int_var( rtn_int_var rtn);

typedef struct {
	int n_conf;
	double mean;
	double *a;
	double error;
	} cluster_jk;

	extern void init_cluster_jk (cluster_jk * c,int n_cluster , int vector_length);
	
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <GL/freeglut.h>
#include "color.h"
#include <gsl/gsl_const_mksa.h>
#include "varie.h"

#define N_SERIES 25
#define R_MIN 0
#define R_MAX 70
#define V_MAX -5
//#define H_BAR (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)
#define NORM  0.1
#define PI 3.14159
#define SIGMA 20

#define D_T 0.5
#define R_PDF 30
/** 
 * width and heigth of the matrix 
 * ad ogni elemento della matrice corrisponde un pixel
 */
const int W = 200;
const int H = 200;
/* reticular pass */
double a = 1;
/* matrix for the wave function */
gsl_matrix_complex* psi ;	
/* for displaying */
int isActive, modeView;
size_t time;
/* Physical parameters*/
double kinetic_constant; 
double harmonic_constant;
gsl_matrix_complex *temp;
//gsl_matrix_complex *matrix_sum;
gsl_matrix_complex *step;
//gsl_matrix_complex *sub_psi;
//gsl_matrix_complex *ham_psi;
/*	
double factorial( int n){
	int i ;
	double result = 1;
	for( i = 1 ; i < n+1; i++)
		result *=  (double) i;
	return ( result);
	}
*/
/* Calcola la norma quadra di Psi. Nota che è ~ \int | \psi |^2, non il prodotto di matrici*/
double matrix_complex_norm ( gsl_matrix_complex *input){
	int w = (int) input->size1;
	int h = (int) input ->size2;
	int i,j;
	double norm_squared = 0.0;
	for( i = 0 ; i < w; i++){
		for(j = 0; j<h; j++){
			norm_squared += 
			pow(gsl_complex_abs(gsl_matrix_complex_get(input,i,j)),2);
		}
	}
	return (sqrt(norm_squared));
}
/* Ho aggiunto un E^ikx */
gsl_complex circular_step_pdf( double x , double y ){
		return gsl_complex_rect ( 40*sin(-1e2*x)*gaussPdf(SIGMA*a,50,x)*gaussPdf(SIGMA*a,0,y),40*cos(-1e2*x)*gaussPdf(SIGMA*a,50,x)*
		gaussPdf(SIGMA*a,0,y));
	}
/*
double V_step_tunnel (double x , double y ){
	}
*/
double potential( int i ,int j  ){
	//return 0.0;
	//return( (W/2/PI)*(W/2/PI)*sin( 2*PI/W*i)*sin( 2*PI/W*i)+(H/2/PI)*(H/2/PI)*cos( 2*PI/H*j)*cos( 2*PI/H*j));
	//return ( 1e-3*harmonic_constant*a*a*(i*i+j*j));
	//return (10);
	
	if ( i < R_MAX && j < R_MAX)
		return -V_MAX;
	else
		return 0.0 ;
	
}

void init_wave_function (gsl_matrix_complex *input , gsl_complex (*pdf) ( double x , double y ) ) {
	int i , j ;
	int w = (int) input->size1;
	int h = (int) input ->size2;
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
		    gsl_matrix_complex_set(input,i,j, pdf(a*(i-w/2),a*(j-h/2) ) );
		}
	}
	gsl_matrix_complex_scale( input,gsl_complex_rect(1.0/matrix_complex_norm(input),0));
	}

/* Set matrix out equal to matrix in */
	void matrix_equal(gsl_matrix_complex *in , gsl_matrix_complex * out){
	int i,j ;
	int w = (int) in->size1;
	int h = (int) in->size2;
	int w1 = (int) out->size1;
	int h1 = (int) out->size2;
	if( w == w1 && h==h1){
	for( i = 0 ; i < h ; i++)
		for( j = 0 ; j < w ; j++)
			gsl_matrix_complex_set(out, i, j, gsl_matrix_complex_get(in, i, j));
	}
	else{
		printf("Matrix_equal : matrices have different dimensions!\n");
		exit(EXIT_FAILURE);
	}
	}
/* Per una |psi> data (in) restituisce H | psi > e la salva in out*/
void  hamiltonian ( gsl_matrix_complex* in , gsl_matrix_complex * out){
	int w = (int) in->size1;
	int h = (int) in ->size2;
	int i, j ;
	gsl_complex sum ;
	for ( i = 0 ; i< h ; i++){
		for( j = 0 ; j< w ; j++){
			sum = GSL_COMPLEX_ZERO;
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,i ,(j+1)%w));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,i ,(j-1 + w)% w));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,(i+1)%h ,j ));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,(i-1 +h)%h ,j ));
			sum = gsl_complex_mul_real(sum , 0.25);
			sum = gsl_complex_add( sum , gsl_matrix_complex_get(in , i , j) );
			gsl_matrix_complex_set(out, i, j,
			  gsl_complex_add(
			    gsl_complex_mul_real(sum,kinetic_constant), 
			    gsl_complex_mul_real(gsl_matrix_complex_get(in,i,j),potential((i)%h -h/2,(j)%w-w/2) ) ) );
		}
	}
		//gsl_matrix_complex_scale(out, gsl_complex_rect(50.0/matrix_complex_norm(in),0));
}
	
	/**
	 * H = p^2 / 2m + V(x)
	 * U(t,0) = exp(- i H t) ~ Id -iH t| psi> + ... -(it H)^n /(n!).... nmax= N_SERIES 
	 */
void compute ( gsl_matrix_complex *input ){
	int i;
	/* Matrix_equal( in, out) */
	//matrix_equal( input, matrix_sum);
	matrix_equal(input, temp);
	  for (i = 1; i < N_SERIES ; i++ ){
	/* hamiltonian calcola H |temp > e la salva in step
		hamiltonian ( in, out)*/
		hamiltonian(temp,step);
	/* Moltiplica la matrice per -i(dt)/n, in modo da ricostruire il fattoriale */
		gsl_matrix_complex_scale(step,gsl_complex_div_real(gsl_complex_rect(0,-D_T), (double) i));
	/* salva step in temp */
		matrix_equal(step ,temp);
		gsl_matrix_complex_add( input, step);
		}
	/* Calcola se la differenza di psi è uguale all'hamiltoniana */
	/*matrix_equal(input,sub_psi);
	gsl_matrix_complex_sub( sub_psi,matrix_sum);
	gsl_matrix_complex_scale(sub_psi, gsl_complex_rect(1.0/ (double) D_T,0));
	hamiltonian(input, ham_psi);
	gsl_matrix_complex_add(sub_psi, ham_psi);
	//printf("%e\n", matrix_complex_norm(sub_psi));
	matrix_equal( matrix_sum, input);
	printf("%e \n", matrix_complex_norm(input));
	*/
	//matrix_equal(matrix_sum,input);
	//gsl_matrix_complex_scale(input, gsl_complex_rect(1.0/matrix_complex_norm(input),0));
	//printf("%e \n", matrix_complex_norm(input));
	}

 
 /** ****************************************
  * 
  * OpenGl stuff
  *******************************************/
void GLInit()
{
    glDisable(GL_DEPTH_TEST);
    glClearColor(0.0 ,0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glOrtho(0,W,0,H,0,1);
}

void displayF()
{
	int w = (int) (psi->size1);
	int h = (int) (psi->size2);
	float* data = (float*)malloc(3*(w*h)*sizeof(float));
	gsl_complex tmp ;
    if(data == NULL)
        printf("\n\nmalloc fail!!!\n\n");
    rgb_t color;
    int k;
    for(k = 0; k < w*h; k++){
		tmp = gsl_matrix_complex_get(psi, k/w,  k%w);
        if(modeView == 1)
            color = d2rgb(gsl_complex_abs(tmp)/ NORM);
        if(modeView == 2)
            color = d2rgb( tmp.dat[0]/NORM) ;
        if(modeView == 3)
            color = d2rgb( tmp.dat[1]/NORM);
        data[3*k+0] = color.r;
        data[3*k+1] = color.g;
        data[3*k+2] = color.b;
	}
	glRasterPos2i(0,0);
    glDrawPixels(W,H,GL_RGB,GL_FLOAT,data);
    glutSwapBuffers();
    free(data);
}

void idleF()
{
    time += isActive;
    if(isActive)
       compute(psi);
    glutPostRedisplay();
}

void keyboardF(unsigned char key, int mouseX, int mouseY)
{
    switch(key)
    {
        case 'q': case 'Q': case 27:
	    gsl_matrix_complex_free( psi );
            exit(EXIT_SUCCESS);
        case ' ':
            isActive =! isActive;
            break;
        case '+':
            compute(psi);
            break;
        case '1':
            modeView = 1;
            break;
        case '2':
            modeView = 2;
            break;
        case '3':
            modeView = 3;
            break;
    }
}
int main (int argc, char *argv[]){
	kinetic_constant = 2;
	harmonic_constant = 1;
	temp = gsl_matrix_complex_alloc (W,H) ;
	//matrix_sum = gsl_matrix_complex_alloc(W,H);
	step = gsl_matrix_complex_alloc(W,H);
	//sub_psi = gsl_matrix_complex_alloc(W,H);
	//ham_psi = gsl_matrix_complex_alloc(W,H);
	psi = gsl_matrix_complex_alloc(W,H);
	init_wave_function( psi , circular_step_pdf );
	time = modeView = isActive = 1;
	//printf("%e \n", factorial(4.0));
		glutInit(&argc, argv);
    glutInitWindowSize((int)W,(int)H);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
    glutCreateWindow("Schrodinger"); 
    GLInit();
    glutDisplayFunc(displayF); 
    glutIdleFunc(idleF);
    //glutMouseFunc(mouseF);
    glutKeyboardFunc(keyboardF);
    glutMainLoop();
	gsl_matrix_complex_free(psi);
	gsl_matrix_complex_free(temp);
	//gsl_matrix_complex_free(matrix_sum);
	gsl_matrix_complex_free(step);
	//gsl_matrix_complex_free(sub_psi);
	//gsl_matrix_complex_free(ham_psi);
		return(EXIT_SUCCESS);
	}

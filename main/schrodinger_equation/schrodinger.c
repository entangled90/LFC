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

#define N_SERIES 30
#define R_MIN 0
#define R_MAX 10
#define V_MAX -1
//#define H_BAR (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)
#define NORM  0.1
#define PI 3.14159
#define SIGMA 5
#define RETICOLO 10

#define D_T 0.7
/** 
 * width and heigth of the matrix 
 * ad ogni elemento della matrice corrisponde un pixel
 */
const int W = 150;
const int H = 150;
/* reticular pass */

/**
 *NOTA BENE
 * è necessario modificare i parametri dell'algoritmo (N_SERIES e D_T) a seconda del potenziale scelto.
 * Per potenziali che danno un contributo piccolo agli autovalori dell'hamiltoniana si può porre N_SERIES ~25 e D_T 0.5
 * Nel caso di potenziali più grandi (in questo caso essenzialmente quelli che valgono ~10/100 (es i muri,diffrazione ecc)
 * è necessario porre N_SERIES ~30/40 e D_T= 0.1
 * In ogni caso è possibile controllare l'unitarietà dell'operatore di evoluzione temporale facendo stampare la norma della matrice
 * decommentando l'istruzione
 * printf("%e \n", matrix_complex_norm(input) - 1);
 * all'ultima riga della funzione "compute()"
 */
/* matrix for the wave function */
gsl_matrix_complex* psi ;	
/* for displaying */
int isActive, modeView;
size_t time;
/* Physical parameters*/
const double kinetic_constant = 2; 
const double harmonic_constant = 1e-4;
gsl_matrix_complex *temp;
gsl_matrix_complex *step;

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

gsl_complex circular_step_pdf( double x , double y ){
		  /*Momento verso il basso*/
  return gsl_complex_rect ( 40*sin(-1e2*x)*gaussPdf(SIGMA,20,x)*gaussPdf(SIGMA,0,y),40*cos(-1e2*x)*gaussPdf(SIGMA,20,x)*
		gaussPdf(SIGMA,0,y));
  /*Gaussiana nell'origine*/
  //return gsl_complex_rect ( gaussPdf(SIGMA,0,x)*gaussPdf(SIGMA,0,y),gaussPdf(SIGMA,0,x)*gaussPdf(SIGMA,0,y));
  
}


double potential( int i ,int j  ){
  return 0.0;
 //  return ( harmonic_constant*(i*i+j*j));
/** Buca sferica*/
/*	
	if ( i*i + j*j < R_MAX*R_MAX )
		return V_MAX;
	else
		return 0.0 ;
*/	
/**Barriera di potenziale orizzontale*/
 /*
  * if( abs(i) < R_MAX)
    return V_MAX;
  else
    return 0.0;
*/
/**Punti ~reticolo*/
/*
if ( (( abs(i) == W/8) && (j == 0)))
    return 1e2;
else if ( (( abs(j) == H/8) && (i == 0)))
    return 1e2;
else if ( i ==0 && j == 0)
    return 1e2;
else
  return 0;
*/
/** Diffrazione doppia fenditura*/
/*
  if( i == 0){
    if (abs(j) < 5)
      return 1e2;
    else if(abs(j) > 10)
      return 1e2;
    else
      return 0;
  }
  else
    return 0;
*/
/** ~ Potenziale periodico */
//return ( 1e-1*(sin (2*PI*RETICOLO*i/(double) W) + sin( 2*PI*RETICOLO*j/(double) H)) );
/** ~ Potenziale periodico  con oscillatore ~*/
//return ( 1e-2*( abs(i)*sin (2*PI*RETICOLO*i/(double) W) + abs(j)*sin( 2*PI*RETICOLO*j/(double) H)) );
/** Reticolo di atomi ~*/
/*
 * if ( ( ( i%RETICOLO) == 0) && ( (j%RETICOLO== 0))){
    return (1e1);
  }
else
  return (1e-1*(sin (2*PI*RETICOLO*i/(double) W) + sin( 2*PI*RETICOLO*j/(double) H)) );
*/
  
}


void init_wave_function (gsl_matrix_complex *input , gsl_complex (*pdf) ( double x , double y ) ) {
	int i , j ;
	int w = (int) input->size1;
	int h = (int) input ->size2;
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
		    gsl_matrix_complex_set(input,i,j, pdf((i-w/2),(j-h/2) ) );
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
}
	
	/**
	 * H = p^2 / 2m + V(x)
	 * U(t,0) = exp(- i H t) ~ Id -iH t| psi> + ... -(it H)^n /(n!).... nmax= N_SERIES 
	 */
void compute ( gsl_matrix_complex *input ){
	int i;
	/* Matrix_equal( in, out) */
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
	
	//gsl_matrix_complex_scale(input, gsl_complex_rect(1.0/matrix_complex_norm(input),0));
	/*Stampa la differenza fra la norma e l'unità→ per controllare l'unitarietà dell'operatore*/
	//printf("%e \n", matrix_complex_norm(input) - 1);
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
            color = d2rgb( GSL_REAL(tmp)*NORM) ;
        if(modeView == 3)
            color = d2rgb( GSL_IMAG(tmp)*NORM);
	if(modeView == 4)
	    color = d2rgb (gsl_complex_arg(tmp)*NORM);
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
	case '4':
	    modeView = 4;
            break;
    }
}
int main (int argc, char *argv[]){
     temp = gsl_matrix_complex_alloc (W,H) ;
    step = gsl_matrix_complex_alloc(W,H);
    psi = gsl_matrix_complex_alloc(W,H);
    init_wave_function( psi , circular_step_pdf );
    time = modeView = isActive = 1;
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

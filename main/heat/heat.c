#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>
#include <math.h>
#include <GL/freeglut.h>
#include "color.h"
#include <gsl/gsl_const_mksa.h>
#include "varie.h"

#define N_SERIES 2
#define X_MAX 50
#define Y_MAX 50
#define NORM 700
#define PI 3.14159
#define SIGMA 1500
#define D_T 1
#define TEMP_EXT 100
#define DELTA_TEMP 0
/** 
 * width and heigth of the matrix 
 * ad ogni elemento della matrice corrisponde un pixel
 */
const int W = 400;
const int H = 400;
/* reticular pass */

/* matrix for the wave function */
gsl_matrix* psi ;	
/* for displaying */
int isActive,modeView;
size_t time;
/* Physical parameters*/
double diffusive_constant; 
gsl_matrix *temp;
gsl_matrix *step;

double circular_step_pdf( double x , double y ){
  return  (5e-2*fabs(x*x*x)*y*y*gaussPdf(SIGMA,0,x)*gaussPdf(SIGMA,0,y));
}

double rectangular_step_pdf(double x, double y){
  if( (fabs(x) < X_MAX)){
    if( (fabs(y) < Y_MAX))
      return 500;
    else
      return 0;
  }
    else
    return 0;
}



void init_wave_function (gsl_matrix *input , double (*pdf) ( double x , double y ) ) {
	int i , j ;
	int w = (int) input->size1;
	int h = (int) input ->size2;
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
		    gsl_matrix_set(input,i,j, pdf((i-w/2),(j-h/2) ) );
		}
	}
	}

/* Set matrix out equal to matrix in */
	void matrix_equal(gsl_matrix *in , gsl_matrix * out){
	int i,j ;
	int w = (int) in->size1;
	int h = (int) in->size2;
	for( i = 0 ; i < h ; i++)
		for( j = 0 ; j < w ; j++)
			gsl_matrix_set(out, i, j, gsl_matrix_get(in, i, j));
	}
/* Per una |psi> data (in) restituisce H | psi > e la salva in out*/
void  hamiltonian ( gsl_matrix* in , gsl_matrix * out){
	int w = (int) in->size1;
	int h = (int) in ->size2;
	int i, j ;
	double  sum ;
	for ( i = 0 ; i< h ; i++){
		for( j = 0 ; j< w ; j++){
			sum = 0.0;
			if((j+1)< w)
			  sum -=gsl_matrix_get(in,i ,(j+1)%w);
			else
			  sum -=TEMP_EXT + DELTA_TEMP;
			if( (j-1)>0)
			sum -= gsl_matrix_get(in,i ,(j-1 + w)% w);
			else
			  sum -=TEMP_EXT+DELTA_TEMP;
			if( (i+1)<h)
			sum -= gsl_matrix_get(in,(i+1)%h ,j);
			else
			  sum -=TEMP_EXT;
			if ((i-1)>0)
			  sum -= gsl_matrix_get(in,(i-1 +h)%h ,j );
			else
			  sum -= TEMP_EXT;
			sum *= 0.25;
			sum += gsl_matrix_get(in , i , j);
			gsl_matrix_set(out, i, j,sum*diffusive_constant); 
		}
	}
}
	
	/**
	 * H = p^2 / 2m + V(x)
	 * U(t,0) = exp(- i H t) ~ Id -iH t| psi> + ... -(it H)^n /(n!).... nmax= N_SERIES 
	 */
void compute ( gsl_matrix *input ){
	int i;
	/* Matrix_equal( in, out) */
	matrix_equal(input, temp);
	  for (i = 1; i < N_SERIES ; i++ ){
	/* hamiltonian calcola H |temp > e la salva in step
		hamiltonian ( in, out)*/
		hamiltonian(temp,step);
	/* Moltiplica la matrice per -(dt)/n, in modo da ricostruire il fattoriale */
		gsl_matrix_scale(step,(-D_T/(double) i));
	/* salva step in temp */
		matrix_equal(step ,temp);
		gsl_matrix_add( input, step);
	  }
	
	//gsl_matrix_scale(input, gsl_complex_rect(1.0/matrix_complex_norm(input),0));
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
	double tmp ;
    if(data == NULL)
        printf("\n\nmalloc fail!!!\n\n");
    rgb_t color;
    int k;
    for(k = 0; k < w*h; k++){
	tmp = gsl_matrix_get(psi, k/w,  k%w);
        color = d2rgb(fabs(tmp)/ NORM);
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
	    gsl_matrix_free( psi );
            exit(EXIT_SUCCESS);
        case ' ':
            isActive =! isActive;
            break;
        case '+':
            compute(psi);
            break;
    }
}
int main (int argc, char *argv[]){
    diffusive_constant = 1;
    temp = gsl_matrix_alloc (W,H) ;
    step = gsl_matrix_alloc(W,H);
    psi = gsl_matrix_alloc(W,H);
    init_wave_function( psi ,rectangular_step_pdf );
    time = modeView = isActive = 1;
    glutInit(&argc, argv);
    glutInitWindowSize((int)W,(int)H);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
    glutCreateWindow("Heat"); 
    GLInit();
    glutDisplayFunc(displayF); 
    glutIdleFunc(idleF);
    glutKeyboardFunc(keyboardF);
    glutMainLoop();
    gsl_matrix_free(psi);
    gsl_matrix_free(temp);
    gsl_matrix_free(step);
    return(EXIT_SUCCESS);
}

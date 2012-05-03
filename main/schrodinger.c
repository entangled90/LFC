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


#define R_MIN 80
#define R_MAX 100
#define V_MAX 10
#define H_BAR (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)
#define D_T 1e-5
#define R_PDF 30
/* width and heigth of the matrix */
const int W = 600;
const int H = 600;
/* reticular pass */
double a = 1e-3;
/* matrix for the wave function */
gsl_matrix_complex * psi ;
/* for displaying */
int isActive, modeView;
size_t time;
/* Physical parameters*/
double m = 1e-28;
double omega = 1e-4;
double kinetic_constant; 
double harmonic_constant; 

gsl_complex circular_step_pdf( double x , double y ){
		return gsl_complex_rect ( gaussPdf(a,200*a,x)*gaussPdf(a,0,y)/100, gaussPdf(a,200*a,x)*gaussPdf(a,0,y)/100);
	}

double V_step_tunnel (double x , double y ){
	if ( ( x*x + y*y > R_MIN*R_MIN*a*a) && (x*x +y*y < R_MAX*R_MAX*a*a))
		return V_MAX;
	else
		return 0 ;
	}

double V_parabolic ( double x , double y ){
	
	return ( harmonic_constant*(x*x+y*y));
	}

void init_wave_function (gsl_matrix_complex *input , gsl_complex (*pdf) ( double x , double y ) ) {
	int i , j ;
	int w = (int) input->size1;
	int h = (int) input ->size2;
	
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
			gsl_matrix_complex_set(input,i,j, pdf( a*(i-w/2),a*(j-h/2) ) );
		}
	}
	}
	
void compute ( gsl_matrix_complex *input ){
	int i, j ,status = 0;
	int w = (int) input->size1;
	int h = (int) input->size2;
	gsl_complex *increment = malloc(sizeof(gsl_complex));
	gsl_complex *tmp = malloc(sizeof(gsl_complex));
	gsl_matrix_complex *increment_matrix = gsl_matrix_complex_alloc (w,h) ;
	gsl_matrix_complex *psi_old = gsl_matrix_complex_alloc (w,h);
	gsl_matrix_complex_memcpy(psi_old, input);
	for ( i = 0 ; i< w ; i++){
		for( j = 0 ; j< h ; j++){
			*increment = GSL_COMPLEX_ZERO;
			/* Mi sa che è il laplaciano che fa cagare! */ 	
			if( (i != 0) && (i != w-1) && (j != 0) && (j != h-1)){
				if( i != 0)
					*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i-1,j)); 
				if( i != w-1)
					*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i+1,j));
				if ( j != 0)
					*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i,j-1));
				if ( j != h-1) 
					*increment = gsl_complex_add ( * increment , gsl_matrix_complex_get( psi_old, i,j+1));
				(*tmp) = gsl_matrix_complex_get( psi_old, i,j);
				*increment = gsl_complex_sub( * increment ,gsl_complex_mul_real( *tmp , 4.0)  ); 
				*increment = gsl_complex_mul_real( *increment, kinetic_constant ) ;
				/* Ora in increment ci sta salvato:
				 * laplaciano*h_bar/2m (c'è 1 h_bar solo perchè ci divido davanti*/
				//*tmp = gsl_complex_mul_real(gsl_matrix_complex_get(psi_old,i,j),(V_step_tunnel(a*(i-w/2),a*(j-h/2))));
				//*increment = gsl_complex_sub (*increment ,*tmp);
				*increment = gsl_complex_mul (*increment,(gsl_complex_rect(0,1)));
			}
			gsl_matrix_complex_set( increment_matrix,i,j,*increment);
		}
	}
	/* the result is stored in input, increment_matrix is unchanged
	 * → if status is non-zero, an error occurred
	 */
	status += gsl_matrix_complex_add( input , increment_matrix);
	if( status != 0 ){
		printf("An error occurred in compute(): gsl_matrix_complex_add() \n ");
		fflush(stdout);
	}
	gsl_matrix_complex_free(psi_old);
	gsl_matrix_complex_free(increment_matrix);
	free(increment);
	free(tmp);
	
	}


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
            color = d2rgb(gsl_complex_abs(tmp));
        if(modeView == 2)
            color = d2rgb( tmp.dat[0]) ;
        if(modeView == 3)
            color = d2rgb( tmp.dat[1]);
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
           	gsl_matrix_complex_free(psi);
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
	kinetic_constant = 1e-3;
	harmonic_constant = 5;
	psi = gsl_matrix_complex_alloc((int)W,(int)H);
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
	return(EXIT_SUCCESS);
	}

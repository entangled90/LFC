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

#define N_SERIES 20
#define R_MIN 80
#define R_MAX 100
#define V_MAX 10
#define H_BAR (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)

#define D_T 1e-5
#define R_PDF 30
/* width and heigth of the matrix */
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

int factorial( int n){
	int i ;
	int result = 1;
	for( i = 1 ; i < n+1; i++)
		result *= i;
	return ( result);
	}

gsl_complex circular_step_pdf( double x , double y ){
		return gsl_complex_rect ( 1000*gaussPdf(100*a,0,x)*gaussPdf(50*a,0,y), 0);
	}

double V_step_tunnel (double x , double y ){
	if ( ( x*x + y*y > R_MIN*R_MIN*a*a) && (x*x +y*y < R_MAX*R_MAX*a*a))
		return V_MAX;
	else
		return 0 ;
	}

double potential( int i ,int j  ){
	
	return ( harmonic_constant*a*a*(i*i+j*j));
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

	
void  hamiltonian ( gsl_matrix_complex* in , gsl_matrix_complex * out){
	int w = (int) in->size1;
	int h = (int) in ->size2;
	int i, j ;
	gsl_complex sum ;
	for ( i = 0 ; i< w ; i++){
		for( j = 0 ; j< h ; j++){
			sum = GSL_COMPLEX_ZERO;
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,i ,(j+1)%w));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,i ,(j-1 + w)% w));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,(i+1)%h ,j ));
			sum = gsl_complex_sub( sum, gsl_matrix_complex_get(in,(i-1 +h)%h ,j ));
			sum = gsl_complex_mul_real(sum , 0.25);
			sum = gsl_complex_sub( sum , gsl_matrix_complex_get(in , i , j) );
			gsl_matrix_complex_set(out, i, j, gsl_complex_add_real(gsl_complex_mul_real(sum,kinetic_constant), potential(i,j)));
		}
	}
	}
	/* Set matrix out equal to matrix in */
	void matrix_equal(gsl_matrix_complex *in , gsl_matrix_complex * out){
	int i,j ;
	int w = (int) in->size1;
	int h = (int) in->size2;
	for( i = 0 ; i < w ; i++)
		for( j = 0 ; j < h ; j++)
			gsl_matrix_complex_set(out, i, j, gsl_matrix_complex_get(in, i, j));

	}
	
void compute ( gsl_matrix_complex *input ){
	int i;
	int w = (int) input->size1;
	int h = (int) input->size2;
	gsl_matrix_complex *temp = gsl_matrix_complex_alloc (w,h) ;
	gsl_matrix_complex *matrix_sum = gsl_matrix_complex_alloc(w,h);
	gsl_matrix_complex *step = gsl_matrix_complex_alloc(w,h);
	
		matrix_equal( input, matrix_sum);
		for (i = 1; i < N_SERIES ; i++ ){
			hamiltonian(temp,step);
			matrix_equal(temp,step);
			gsl_matrix_complex_scale(step,gsl_complex_mul_real(gsl_complex_pow_real(gsl_complex_rect(0,D_T),i),factorial(i)));
			gsl_matrix_complex_add( matrix_sum, step);
			}
		gsl_matrix_complex_add( input, matrix_sum);
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
	kinetic_constant = 1e-2;
	harmonic_constant = 100;
	psi = gsl_matrix_complex_alloc(W,H);
	init_wave_function( psi , circular_step_pdf );
	time = modeView = isActive = 1;
	printf("%d \n", factorial(4));
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

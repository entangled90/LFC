#include <stdio.h>
#include <stdlib.h>
#include <types.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <GL/freeglut.h>
#include "color.h"
#include <gsl/gsl_const_mksa.h>


#define R_MIN 0.4
#define R_MAX 0.5
#define V_MAX 0.2
#define H_BAR GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
#define D_T 0.0001
#define R_PDF 0.3
/* width and heigth of the matrix */
int w = 500;
int h = 500;
/* reticular pass */
double a = 1;
/* matrix for the wave function */
gsl_matrix_complex * psi ;
/* for displaying */
int isActive, modeView;
size_t time;



gsl_complex circular_step_pdf( double x , double y ){
	if ( ( x*x + y*y < R_PDF))
		return gsl_complex_rect ( sqrt(M_PI * R_PDF* R_PDF), 0);
	else
		return gsl_complex_rect(0,0);
	}

double V_step_tunnel (double x , double y ){
	if ( ( x*x + y*y > R_MIN*R_MIN ) && (x*x +y*y < R_MAX*R_MAX))
		return V_MAX;
	else
		return 0 ;
	}

void init_wave_function (gsl_matrix_complex *psi , gsl_complex (*pdf) ( double x , double y ) ) {
	int i , j ;
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
			gsl_matrix_complex_set(psi,i,j, pdf( a*i  , a*j ) );
		}
	}
	}
	
void compute ( gsl_matrix_complex *input ){
	int i, j ,status = 0;
	gsl_complex c1,c2;
	gsl_complex *increment = &(c1);
	gsl_complex *tmp = &(c2);
	gsl_matrix_complex *increment_matrix = gsl_matrix_complex_alloc ( w,h) ;
	gsl_matrix_complex *psi_old = gsl_matrix_complex_alloc (w,h);
	gsl_matrix_complex_memcpy(psi_old, input);
	for ( i = 1 ; i<= w ; i++){
		for( j = 1 ; j<= h ; j++){
			*increment = GSL_COMPLEX_ZERO;
			if( i > 1)
				*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i-1,j)); 
			if( i < w)
				*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i+1,j));
			if ( j > 1)
				*increment = gsl_complex_add ( *increment , gsl_matrix_complex_get( psi_old, i,j-1));
			if ( j < h) 
				*increment = gsl_complex_add ( * increment , gsl_matrix_complex_get( psi_old, i,j+1));
			(*tmp) = gsl_matrix_complex_get( psi_old, i,j);
			*increment = gsl_complex_sub( *increment ,gsl_complex_mul( gsl_complex_rect(4,0),*tmp)); 
			GSL_SET_COMPLEX(tmp, a*a*H_BAR,0);
			*increment = gsl_complex_div( *increment, *tmp);
			GSL_SET_COMPLEX(tmp,V_step_tunnel( a*i, a*j)/((double) H_BAR),0);
			*increment = gsl_complex_sub( *increment , *tmp);
			*increment = gsl_complex_mul ( *increment, gsl_complex_rect(D_T,1)); 
			gsl_matrix_complex_set( increment_matrix,i,j, *increment);
		}
	}
	gsl_matrix_complex_free(psi_old);
	gsl_matrix_complex_free(increment_matrix);
	/* the result is stored in input, increment_matrix is unchanged
	 * → if status is non-zero, an error occurred
	 */
	status += gsl_matrix_complex_add( input , increment_matrix);
	if( status != 0 )
		printf("An error occurred in compute() \n ");
	}


void GLInit()
{
    glDisable(GL_DEPTH_TEST);
    glClearColor(0.0 ,0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glOrtho(0,w,0,h,0,1);
}

void displayF()
{
    float* data = (float*)malloc(3*(w*h)*sizeof(float));
	gsl_complex tmp ;
    if(data == NULL)
        printf("\n\nmalloc fail!!!\n\n");
    rgb_t color;
    int k;
    for(k = 0; k < w*h; k++){
		tmp = gsl_matrix_complex_get(psi, k/w,  k%w);
        if(modeView == 1)
            color = d2rgb(gsl_complex_abs( tmp ));
        if(modeView == 2)
            color = d2rgb( tmp.dat[0]) ;
        if(modeView == 3)
            color = d2rgb( tmp.dat[1]);
        data[3*k+0] = color.r;
        data[3*k+1] = color.g;
        data[3*k+2] = color.b;
    }
    glRasterPos2i(0,0);
    glDrawPixels(w,h,GL_RGB,GL_FLOAT,data);
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
	psi = gsl_matrix_complex_alloc (w,h);
	init_wave_function( psi , circular_step_pdf );
	time = modeView = isActive = 1;
    glutInit(&argc, argv);
    glutInitWindowSize(w,h); 
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

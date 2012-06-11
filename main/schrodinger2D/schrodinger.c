#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <GL/glut.h>
#include "color.h"
#include "varie.h"
#include "schrodinger_lib.h"

#define NORM  0.1
#define PI 3.14159
/** 
 * width and heigth of the matrix 
 * ad ogni elemento della matrice corrisponde un pixel
 */
const int W = 150;
const int H = 150;
/* reticular pass */

/* matrix for the wave function */
gsl_matrix_complex* psi ;	
/* for displaying */
int isActive, modeView;
size_t time;
/* Physical parameters*/
gsl_matrix_complex *temp;
gsl_matrix_complex *step;


 
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
       compute(psi,temp,step);
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
            compute(psi,temp,step);
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

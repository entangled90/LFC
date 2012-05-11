#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "color.h"
#include <GL/freeglut.h>
#define H 0.001
/* Friction coefficient */
#define Q 0.1
/* Extern force amplitude */
#define B -0.5
#define OMEGA_EXT (2.0/3.0)
#define N_STEPS 100000000
#define PI 3.14159265
int N_POINTS_DISPLAYED= 1e4;
int n_points = 1;
float a =1 ;
float ds;
float scalef;
int ACTIVE;

double k1[4];
double k2[4];
double t;

 struct point2D  {
  double x;
  double y;
	struct point2D *next;
} ;
struct point2D b ;
struct point2D * head ;
struct point2D * tail;
/*
point2D *seek_element ( point2D *previous){
	if ( previos == tail)
		return ;
	return ( seek_element(previos->next));
	}
*/
struct point2D * get_element ( struct point2D head, int position  ) {
	int i = 0;
	struct point2D *cursor;
	cursor = &head;
	for( i = 0; i< position; i++){
		cursor = cursor->next;
	}
	return ( cursor);
	}
void free_list ( struct point2D *head , struct point2D *tail, int n){
	struct point2D cursor, *tmp;
	int i ;
	cursor = *head;
	tmp = head;
	for ( i = 0 ; i<n ; i++){
		cursor = *(tmp);
		free(tmp);
		tmp=cursor.next;
		
	}
	}

/* Algoritmo per runge-kutta in 2D */
void kn_fill (double *k1_vec, double *k2_vec, struct point2D p , double t , double (*f1) (double ,double,double ), double (*f2) (double, double,double)){
  int i = 0;
  k1_vec[0]= H*f1(p.x ,p.y ,t);
  k2_vec[0]= H*f2(p.x, p.y ,t);
  for(i = 0; i< 2 ; i++){
    k1_vec[i+1]= H*f1(p.x+k1_vec[i]/2,p.y,t+H/2 ); 
    k2_vec[i+1]= H*f2(p.x,p.y+k2_vec[i]/2,t+H/2 ); 
  }
  k1_vec[3]= H * f1(p.x+k1_vec[2],p.y,t+H);
  k2_vec[3]= H * f2(p.x,p.y+k1_vec[2],t+H);

}

double f2 ( double x1, double x2, double time){
  return ( 4*(1-x1*x1)*x2-x1);
  //-sin(x1) -Q*x2+B*cos(OMEGA_EXT*time));
}
double f1 ( double x1, double x2, double time){
  return (x2);

}
double increment_k ( double * k ){
  return ( k[0]/3+k[1]/6+k[2]/6+k[3]/3);
}

void compute ( ){
	tail->next =malloc(sizeof(struct point2D ));
 kn_fill(k1,k2,	*(tail),t+n_points*H,f1,f2);
 (tail->next)->x = tail->x+increment_k(k1);
 (tail->next)->y = tail->y + increment_k(k2);
 n_points++;
 tail = tail->next;
 }

/* OpenGL stuff */

void drawCircle(struct point2D p,double r)
{
	int i;
	rgb_t color;
	color = d2rgb(fabs( tail->x * tail->y));
	glColor3d(color.r,color.g,color.b);
    glBegin(GL_POLYGON);
    for (i = 0; i < 360; i++)
        glVertex2d(p.x+r*cos(i*PI/180.0),p.y+r*sin(i*PI/180.0));
    glEnd();
}


void drawLine( struct point2D *head )
{
  int i = 0;
	struct point2D * cursor;
    glBegin(GL_LINE_STRIP);
	if( n_points < N_POINTS_DISPLAYED)
		cursor = head;
	else
		cursor = get_element( *head, n_points - N_POINTS_DISPLAYED);
	do{
		rgb_t color;
        color = d2rgb(fabs( cursor->y*cursor->y));
        glColor3d(color.r,color.g,color.b);
        glVertex2d((cursor+i)->x,(cursor+i)->y);
		cursor = cursor->next;
    }while ( cursor != tail );
    glEnd();
}

void displayF()
{
    glClear(GL_COLOR_BUFFER_BIT);
    drawLine(head);
    drawCircle(*tail ,0.01);
    glutSwapBuffers();
}

void reshapeF(int w,int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-w*a/h+b.x,w*a/h+b.x,-a+b.y,a+b.y,1,-1);
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0.5,0.5,0.5,0.0);
}

void init()
{
  ACTIVE = 1;
  ds = 0.05;
  a = 1;
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    reshapeF(viewport[2],viewport[3]);
}

void keyboardF(unsigned char key, int x, int y)
{
    switch(key)
    {
        case '+':
            glTranslatef(b.x,b.y,0);
            glScalef(scalef,scalef,1);
            glTranslatef(-b.x,-b.y,0);
            ds/=scalef;
            a/=scalef;
            break;
        case '-':
            glTranslatef(b.x,b.y,0);
            glScalef(1/scalef,1/scalef,1);
            glTranslatef(-b.x,-b.y,0);
            ds*=scalef;
            a*=scalef;
            break;
        case 'p': case 'P': case ' ':
            ACTIVE=!ACTIVE;
            break;
        case 'f': case 'F':
            glutFullScreenToggle();
            break;
        case 'r': case 'R':
            init();
            break;
        case 't':
			if ( N_POINTS_DISPLAYED< pow(2,30))
				N_POINTS_DISPLAYED *=2;
			break;
		case 'e':
			if(N_POINTS_DISPLAYED > 4)
				N_POINTS_DISPLAYED /=2;
			break;
        case 'q': case 'Q': case 27:
            exit(0);
            break;
    }
}

void specialKeyboardF(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_F11:
            glutFullScreenToggle();
            break;
        case GLUT_KEY_UP:
            glTranslatef(0,-ds,0);
            b.y+=ds;
            break;
        case GLUT_KEY_DOWN:
            glTranslatef(0,ds,0);
            b.y-=ds;
            break;
        case GLUT_KEY_LEFT:
            glTranslatef(ds,0,0);
            b.x-=ds;
            break;
        case GLUT_KEY_RIGHT:
            glTranslatef(-ds,0,0);
            b.x+=ds;
            break;
    }
}


void idleF(void)
{
  if(ACTIVE){
      compute();
  }
    glutPostRedisplay();
}


/* End of OpenGL stuff */



int main (int argc, char* argv[]){
  /*  double tmp1,tmp2;
FILE *fp = fopen("../data/differential_equation/runge_kutta.dat","w");
for( i = 0; i< N_STEPS ; i++){
   
     
     tmp1 =increment_k(k1);
    tmp2 = increment_k(k2);
    if( i%10 == 0)
    fprintf(fp,"%e \t %e\n", x1,x2);
  
    } 
      fclose(fp);
  */
  head =malloc (sizeof(struct point2D));
  head->next = head;
  tail = head;
  b.x = 0;
  b.y =0 ;
  t =0;
  scalef=2;
  n_points =1 ;
  head->x= 0.5;
  head->y = 1;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowSize(500, 500);
  glutCreateWindow("Chaotic pendulum");
  init();
  glutDisplayFunc(displayF);
  glutIdleFunc(idleF);
  glutKeyboardFunc(keyboardF);
  glutSpecialFunc(specialKeyboardF);
    //glutMouseFunc(mouseF);
  glutReshapeFunc(reshapeF);
  glutMainLoop();
  free_list(head,tail,n_points);
    return(EXIT_SUCCESS);

}

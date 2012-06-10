#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "varie.h"

//#define PIPE 

/**
 * N DEVE ESSERE MULTIPLO DI 100. altrimenti ci sono problemi con la visualizzazione
 */
#define N 200
#define N_SERIES 30
#define R_MIN 0
#define R_MAX 10
#define V_MAX -1
#define NORM  0.1
#define PI 3.14159
#define VARIANCE 10
#define RETICOLO 10

#define D_T 0.4

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


gsl_matrix_complex* psi ;
double initial_max ;
double z_axis_value;
const double kinetic_constant = 2; 
const double harmonic_constant = 5e-3;
gsl_matrix_complex *temp;
gsl_matrix_complex *step;
int gl_time_prec;
int timediff;
double normalization;
float angle_prec;
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
  return gsl_complex_rect ( 40*sin(-1e2*x)*gaussPdf(VARIANCE,20,x)*gaussPdf(VARIANCE,0,y),40*cos(-1e2*x)*gaussPdf(VARIANCE,20,x)*
		gaussPdf(VARIANCE,0,y));
  
  /*Gaussiana nell'origine*/
//  return gsl_complex_rect ( gaussPdf(VARIANCE,N/12,x)*gaussPdf(VARIANCE,0,y),gaussPdf(VARIANCE,0,x)*gaussPdf(VARIANCE,0,y));
  
}


double potential( int i ,int j  ){
  // return 0.0;
  // return ( -harmonic_constant*(i*i+j*j) );
/** Buca sferica*/

	if ( i*i + j*j < R_MAX*R_MAX )
		return  V_MAX;
	else
		return 0.0 ;

	/**Barriera di potenziale orizzontale*/
	/*
  if( abs(i) < R_MAX)
    return V_MAX;
  else
    return 0.0;
*/
/**Punti ~reticolo*/
/*
if ( (( abs(i) == N/8) && (j == 0)))
    return 1e1;
else if ( (( abs(j) == N/8) && (i == 0)))
    return 1e1;
else if ( i ==0 && j == 0)
    return 1e1;
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
//return ( 1e-1*(sin (2*PI*RETICOLO*i/(double) N) + sin( 2*PI*RETICOLO*j/(double) N)) );
/** ~ Potenziale periodico  con oscillatore ~*/
//return ( 1e0*(-cos(2*PI*RETICOLO*i/(double) N) -cos( 2*PI*RETICOLO*j/(double) N)) );
/** Reticolo di atomi ~*/
/*
 if ( ( ( i%RETICOLO) == 0) && ( (j%RETICOLO== 0))){
    return (1e1);
  }
else
  return (1e-1*(sin (2*PI*RETICOLO*i/(double) N) + sin( 2*PI*RETICOLO*j/(double) N)) );
*/
  
}


double init_wave_function (gsl_matrix_complex *input , gsl_complex (*pdf) ( double x , double y ) ) {
	int i , j ;
	int w = (int) input->size1;
	int h = (int) input ->size2;
	for ( i = 0 ; i < w ; i++){
		for ( j = 0 ; j < h ; j++ ) {
		    gsl_matrix_complex_set(input,i,j, pdf((i-w/2),(j-h/2) ) );
		}
	}
	gsl_matrix_complex_scale( input,gsl_complex_rect(1.0/matrix_complex_norm(input),0));
	// deve ritornare il massimo della funzione d'onda!
	double max = 0;
	for(int i = 0; i < N; i++) {
	  for(int j = 0; j < N; j++) {
	    double z = 0;
	    z = gsl_complex_abs(gsl_matrix_complex_get(psi,i,j));
	    if(fabs(z)>max)
	      max = fabs(z);
	   }
	}
	return max;

}

/* Set matrix out equal to matrix in */
	void matrix_equal(gsl_matrix_complex *in , gsl_matrix_complex * out){
	int i,j ;
	int w = (int) in->size1;
	int h = (int) in->size2;
	for( i = 0 ; i < h ; i++){
		for( j = 0 ; j < w ; j++)
			gsl_matrix_complex_set(out, i, j, gsl_matrix_complex_get(in, i, j));
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
/**
 *  ROBA OPENGL per il 3D!
 */

GLuint program;
GLint attribute_coord2d;
GLint uniform_vertex_transform;
GLint uniform_texture_transform;
GLint uniform_color;
GLuint texture_id;
GLint uniform_mytexture;

float offset_x = 0.0;
float offset_y = 0.0;
float scale = 1.0;

bool interpolate = true;
bool clamp = false;
bool rotate = false;
bool polygonoffset = true;
bool active = true;

int modeView = 0;
int width = 640;
int height = 480;

GLuint vbo[3];

struct point {
  GLfloat x;
  GLfloat y;
};

GLbyte graph[N][N];
unsigned char *frame = NULL;

void savePPM(unsigned char *frame)
{
    FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    int i,j;
    for(i = height-1; i >= 0; i--)
        for(j = 0; j < width; j++)
            fwrite(&frame[(i*width+j)*3], sizeof(unsigned char), 3, f);
    fclose(f);
}
void setupTexture()
{
  // Create our datapoints, store it as bytes
  /*
  double max = 0;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      double z = 0;
      switch(modeView){
        case 0:
          z = gsl_complex_abs(gsl_matrix_complex_get(psi,i,j));
          break;
        case 1:
          z = GSL_REAL(gsl_matrix_complex_get(psi,i,j));
          break;
        case 2:
          z = GSL_IMAG(gsl_matrix_complex_get(psi,i,j));
          break;
      }
      if(fabs(z)>max)
        max = fabs(z);
    }
  }
  */
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      double z = 0;
      switch(modeView){
        case 0:
          z = gsl_complex_abs(gsl_matrix_complex_get(psi,i,j));
          break;
        case 1:
          z = GSL_REAL(gsl_matrix_complex_get(psi,i,j));
          break;
        case 2:
          z = GSL_IMAG(gsl_matrix_complex_get(psi,i,j));
          break;
      }
      graph[i][j] = roundf(normalization*z/initial_max*127 + 128); // c'era z/max prima
    }
  }

  /* Upload the texture with our datapoints */
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexImage2D(
     GL_TEXTURE_2D,      // target
     0,                  // level, 0 = base, no minimap,
     GL_LUMINANCE,       // internalformat
     N,                  // width
     N,                  // height
     0,                  // border, always 0 in OpenGL ES
     GL_LUMINANCE,       // format
     GL_UNSIGNED_BYTE,   // type
     graph
  );

  // Create an array for 101 * 101 vertices
  point vertices[N+1][N+1];

  for(int i = 0; i < N+11; i++) {
    for(int j = 0; j < N+1; j++) {
      vertices[i][j].x = (j - N/2) / (double) N;
      vertices[i][j].y = (i - N/2) / (double) N;
    }
  }

  // Tell OpenGL to copy our array to the buffer objects
  glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);

  // Create an array of indices into the vertex array that traces both horizontal and vertical lines
  GLushort indices[N * (N+1) * 6];
  int i = 0;

  for(int y = 0; y < N; y++) {
    for(int x = 0; x < N; x++) {
      indices[i++] = y * (N+1)+ x;
      indices[i++] = y * (N+1) + x + 1;
    }
  }

  for(int x = 0; x < N+1; x++) {
    for(int y = 0; y < N; y++) {
      indices[i++] = y * (N+1)+ x;
      indices[i++] = (y + 1) * (N+1)+ x;
    }
  }

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[1]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, N * (N+1) * 4 * sizeof *indices, indices, GL_STATIC_DRAW);

  // Create another array of indices that describes all the triangles needed to create a completely filled surface
  i = 0;

  for(int y = 0; y < N+1; y++) {
    for(int x = 0; x < N; x++) {
      indices[i++] = y * (N+1) + x;
      indices[i++] = y * (N+1) + x + 1;
      indices[i++] = (y + 1) * (N+1) + x + 1;

      indices[i++] = y * (N+1) + x;
      indices[i++] = (y + 1) * (N+1) + x + 1;
      indices[i++] = (y + 1) * (N+1) + x;
    }
  }

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[2]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof indices, indices, GL_STATIC_DRAW);

}

int init_resources()
{

  int vertex_texture_units;
  glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &vertex_texture_units);
  if(!vertex_texture_units) {
    fprintf(stderr, "Your graphics cards does not support texture lookups in the vertex shader!\n");
    return 0;
  }

  GLint link_ok = GL_FALSE;

  GLuint vs, fs;
  if ((vs = create_shader("graph.v.glsl", GL_VERTEX_SHADER))   == 0) return 0;
  if ((fs = create_shader("graph.f.glsl", GL_FRAGMENT_SHADER)) == 0) return 0;

  program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    return 0;
  }

  const char* attribute_name;
  attribute_name = "coord2d";
  attribute_coord2d = glGetAttribLocation(program, attribute_name);
  if (attribute_coord2d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }

  const char* uniform_name;
  uniform_name = "vertex_transform";
  uniform_vertex_transform = glGetUniformLocation(program, uniform_name);
  if (uniform_vertex_transform == -1) {
    fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
    return 0;
  }

  uniform_name = "texture_transform";
  uniform_texture_transform = glGetUniformLocation(program, uniform_name);
  if (uniform_texture_transform == -1) {
    fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
    return 0;
  }

  uniform_name = "color";
  uniform_color = glGetUniformLocation(program, uniform_name);
  if (uniform_color == -1) {
    fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
    return 0;
  }

  // Generate texture name
  glActiveTexture(GL_TEXTURE0);
  glGenTextures(1, &texture_id);
  // Create two vertex buffer objects
  glGenBuffers(3, vbo);

  setupTexture();

  return 1;
}

void display()
{
  if(active){
    compute(psi);
    setupTexture();
  }
  
#ifdef PIPE
      frame = (unsigned char*)malloc(3*width*height*sizeof(unsigned char));
      glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,frame);
      fwrite(frame, sizeof(char), 3*(N*N), stdout);
      free(frame);
#endif
  glUseProgram(program);
  glUniform1i(uniform_mytexture, 0);

  glm::mat4 model;

  if(rotate){
    timediff = GLUT_ELAPSED_TIME-gl_time_prec; 
    angle_prec += (float) (timediff*360.0)/ 40000.0f;
    model = glm::rotate(glm::mat4(1.0f), angle_prec, glm::vec3(0.0f, 0.0f, 1.0f));
  //  fprintf(stderr,"%d\n",timediff);
  }
    else
    model = glm::mat4(1.0f);  

  glm::mat4 view = glm::lookAt(glm::vec3(0.0, z_axis_value, 1.8), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 2.0));
  glm::mat4 projection = glm::perspective(45.0f, 1.0f*640/480, 0.1f, 10.0f);

  glm::mat4 vertex_transform = projection * view * model;
  glm::mat4 texture_transform = glm::translate(glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, 1)), glm::vec3(offset_x, offset_y, 0));

  glUniformMatrix4fv(uniform_vertex_transform, 1, GL_FALSE, glm::value_ptr(vertex_transform));
  glUniformMatrix4fv(uniform_texture_transform, 1, GL_FALSE, glm::value_ptr(texture_transform));

  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* Set texture wrapping mode */
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, clamp ? GL_CLAMP_TO_EDGE : GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, clamp ? GL_CLAMP_TO_EDGE : GL_REPEAT);

  /* Set texture interpolation mode */
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, interpolate ? GL_LINEAR : GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, interpolate ? GL_LINEAR : GL_NEAREST);

  /* Draw the triangles, a little dark, with a slight offset in depth. */
  GLfloat grey[4] = {0.5, 0.5, 0.5, 1};
  glUniform4fv(uniform_color, 1, grey);

  glEnable(GL_DEPTH_TEST);

  if(polygonoffset) {
    glPolygonOffset(1, 0);
    glEnable(GL_POLYGON_OFFSET_FILL);
  }

  glEnableVertexAttribArray(attribute_coord2d);
  glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
  glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[2]);
  glDrawElements(GL_TRIANGLES, N * N * 6, GL_UNSIGNED_SHORT, 0);

  glPolygonOffset(0, 0);
  glDisable(GL_POLYGON_OFFSET_FILL);

  /* Draw the grid, very bright */
  GLfloat bright[4] = {2, 2, 2, 1};
  glUniform4fv(uniform_color, 1, bright);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[1]);
  glDrawElements(GL_LINES, N * (N+1) * 4, GL_UNSIGNED_SHORT, 0);

  /* Stop using the vertex buffer object */
  glDisableVertexAttribArray(attribute_coord2d);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glutSwapBuffers();
}

void special(int key, int x, int y)
{
  switch(key) {
    case GLUT_KEY_F1:
      interpolate = !interpolate;
      break;
    case GLUT_KEY_F2:
      clamp = !clamp;
      break;
    case GLUT_KEY_F3:
      rotate = !rotate;
      angle_prec = 0.0;
      break;
    case GLUT_KEY_F4:
      polygonoffset = !polygonoffset;
      break;
    case GLUT_KEY_F5:
      modeView = 0;
      break;
    case GLUT_KEY_F6:
      modeView = 1;
      break;
    case GLUT_KEY_F7:
      modeView = 2;
      break;
    case GLUT_KEY_F8:
      frame = (unsigned char*)malloc(3*width*height*sizeof(float));
      glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,frame);
     savePPM(frame);
     free(frame);
      break;
    case GLUT_KEY_F11:
      glutFullScreenToggle();
      break;
    case GLUT_KEY_LEFT:
      offset_x -= 0.03;
      break;
    case GLUT_KEY_RIGHT:
      offset_x += 0.03;
      break;
    case GLUT_KEY_UP:
      offset_y += 0.03;
      break;
    case GLUT_KEY_DOWN:
      offset_y -= 0.03;
      break;
    case GLUT_KEY_PAGE_UP:
      scale *= 1.5;
      break;
    case GLUT_KEY_PAGE_DOWN:
      scale /= 1.5;
      break;
    case GLUT_KEY_HOME:
      offset_x = 0.0;
      offset_y = 0.0;
      scale = 1.0;
      break;
  }
  glutPostRedisplay();
}

void free_resources()
{
  gsl_matrix_complex_free(psi);
  gsl_matrix_complex_free(temp);
  gsl_matrix_complex_free(step);
  
  glDeleteProgram(program);
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
        case 'q': case 'Q': case 27:
            free_resources();
            exit(EXIT_SUCCESS);
	    break;
        case ' ':
            active =! active;
	    break;
	case '+':
	    normalization*=1.3;
	    break;
	case '-':
	    normalization/=1.3;
	    break;
	case 'w':
	    z_axis_value/=1.2;
	    break;
	case 's':
	  if(z_axis_value>-3.0)
	    z_axis_value*=1.2;
            break;
	case 'i':
	  normalization =1;
	  break;
      
    }
}

int main(int argc, char* argv[]) {
    step = gsl_matrix_complex_alloc(N,N);
  psi = gsl_matrix_complex_alloc(N,N);
  temp = gsl_matrix_complex_alloc(N,N);
  initial_max =init_wave_function( psi , circular_step_pdf );
  normalization = 1;
  z_axis_value = -1.0;
  gl_time_prec = 0;
  timediff=0;
  angle_prec = 0.0f;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA|GLUT_DEPTH|GLUT_DOUBLE);
  glutInitWindowSize(width, height);
  glutCreateWindow("graph");

  GLenum glew_status = glewInit();
  if (GLEW_OK != glew_status) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    return 1;
  }

  if (!GLEW_VERSION_2_0) {
    fprintf(stderr, "No support for OpenGL 2.0 found\n");
    return 1;
  }

  GLint max_units;
  glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &max_units);
  if(max_units < 1) {
	  fprintf(stderr, "Your GPU does not have any vertex texture image units\n");
	  return 1;
  }

  if (init_resources()) {
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutSpecialFunc(special);
    glutKeyboardFunc(keyboard);
    gl_time_prec = glutGet(GLUT_ELAPSED_TIME);
    glutMainLoop();
  }

  free_resources();
  return 0;
} 

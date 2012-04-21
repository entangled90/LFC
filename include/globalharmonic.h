#ifdef MAIN_PROGRAM

#define EXTERN 

#else


#define EXTERN extern
#endif

#define M  1
#define W 1
#define DELTA  2.0
#define Nx 32
#define N_SWEEP 100000000
#define BIN_WIDTH (5*10)
#define THERM_CONST 500
#define N_BIN_CHECK ((N_SWEEP-THERM_CONST)%BIN_WIDTH)
#define K_MAX 5
#define K_START 2
#define N_BIN ((N_SWEEP-THERM_CONST)/BIN_WIDTH)

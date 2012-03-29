 
Routine : dentro modules-> altre directory.

	es:
	static double sum ( double a , double b){
	questa routine è visibile solo nel file in cui è definita, per l'attributo static.
	Non possono essere usate da altri file
}
      le routine visibili dall'esterno NON DEVONO AVERE STATIC davanti es:
	double trapezi ();

      Nei file in cui la si vuole usare, bisogna utilizzare l'include corrispondente.
      nell'include DEVE esserci extern davanti alla definizione della funzione. Vedi random.h

modules/integrazione simpson.c
		     trapezi.c
      

main/ 	run1.c (es)
	

include/	simpson.h	Un include per ogni cartella di moduli
		trapezi.h  
REGOLA:
per ogni directory dei moduli c'è 1 .h. In ogni .h c'è 1 ifndef per ogni file hce c'è dentro i moduli.
dentro gli infdef ci sono extern per le routine visibili in quel file (non static).



MAIN_PROGRAM: Nel programma main, all'inizio definisci "#define MAIN_PROGRAM"


Variabili globali:
		  del tree
		  PERICOLOSISSIME.
		  tutte le variabili che si trovano nelle intestazioni di tutti i file .c per cui non hai messo davanti static.
		  senza static è globale ovunque.
		  Non vengono mai definite in testa ai file, si mettono in global.h in questo modo:
		  
#ifdefined MAIN_PROGRAM

#define EXTERN 

#else

#define EXTERN extern

#endif

EXTERN double pluto;



		  Metti il global dove ti servono le variabili globali.
		  Come funziona la macro:
		  se il file è MAIN_PROGRAM definisce la variabile normalmente : EXTERN = "";
		  se il file è un modulo, EXTERN = extern -> i moduli la vedono come variabile esterna.

		  main_program:
		   static fuori all'inizio del file, le usa solo il main.
		  
		  di un modulo
		    variabile di uso per tutta la routine:
		      in cima al file, prima di definire le routine:
		      static int i;
			  static void pippo (){}
		      NON E' GLOBALE DEL TREE
		      NON DEVI DEFINIRLA SENZA STATIC
		  


		  di una routine

		
		
		
		/* MEMORIA DINAMICA */
		
Memoria statica : + piccola di quella a disposizione in RAM.
  se definisci double rd[1000000000000]; -> non ci sta dà segmentation fault.
  
Memoria assegnata dinamicamente: (da 5e5 in poi)
si dichiara:
double *rd;
/* Malloc accetta una variabile come argomento */
rd = malloc( dim * sizeof(double)); /* nel caso voglio double*/
rd = malloc (dim* sizeof(*rd)); /* + elegante*/
poi lo usi:
PRIMA DI TUTTO METTI TUTTO A ZERO: free non cancella la memoria, la libera all'utilizzo e basta
a=rd[i],
a= *(rd+i);


alla fine della routine:
free(rd); /* Libera la memoria alla ifne della chiamata*/
/*Errore grave dimenticarsi IL FREEEEEEEEEEEEEEE*/









/* Gaussian elimination without pivoting.
Algorithms:
1) At processor 0, complete the data Initiablization including matrix A and B and X
2) Processor0 broadcast data of A, B and X to all other processors 
3) Interleave scheduling, for each processor, it broadcast the specific rows to others processors
4) Wait until all the rows finished calculation, then perform back substitution 
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include "mpi.h"
/*#include <ulocks.h>
#include <task.h>
*/

char *ID;

/* Program Parameters */
#define MAXN 10000  /* Max value of N */
int N;  /* Matrix size */
int procs;  /* Number of processors to use */

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int submit = 0;  /* = 1 if submission parameters should be used */
  int seed = 0;  /* Random seed */
  char uid[L_cuserid + 2]; /*User name */

  /* Read command-line arguments */
  //  if (argc != 3) {
  if ( argc == 1 && !strcmp(argv[1], "submit") ) {
    /* Use submission parameters */
    submit = 1;
    N = 4;
    procs = 2;
    printf("\nSubmission run for \"%s\".\n", cuserid(uid));
      /*uid = ID;*/
    strcpy(uid,ID);
    srand(randm());
  }
  else {
    if (argc == 3) {
      seed = atoi(argv[3]);
      srand(seed);
      printf("Random seed = %i\n", seed);
    }
    else {
      printf("Usage: %s <matrix_dimension> <num_procs> [random seed]\n",
	     argv[0]);
      printf("       %s submit\n", argv[0]);
      exit(0);
    }
  }
    //  }
  /* Interpret command-line args */
  if (!submit) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
    procs = atoi(argv[2]);
    if (procs < 1) {
      printf("Warning: Invalid number of processors = %i.  Using 1.\n", procs);
      procs = 1;
    }
  }

  /* Print parameters */
  printf("\nMatrix dimension N = %i.\n", N);
  printf("Number of processors = %i.\n", procs);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 10) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

int main(int argc, char **argv) {
  /* Timing variables */
  ID = argv[argc-1];
  argc--;
  
  //MPI INIT
  MPI_Init(&argc,&argv);
  MPI_Commm_size(MPI_COMM_WORLD,&procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /* Process program parameters */
  parameters(argc, argv);

  /* Initialize A and B for step 1 of Algorithms*/
  if(myid == 0){
  initialize_inputs();
  /* Print input matrices */
  print_inputs();
  }
  /* Gaussian Elimination */
  gauss();
  MPI_Finalize();
}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */

void gauss() {
  int norm, row, col,i;  /* Normalization row, and zeroing
			* element row and col */
  float multiplier;
  MPI_Status status;
  /* Calculate the Running time */
  double startTime = 0.0;
  double endTime = 0.0;
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid == 0){
    startTime = MPI_Wtime();
  }
 /* Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++) {
    {
      MPI_Bcast(&A[norm][0], N, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&B[norm],1,MPI_FLOAT, 0, MPI_COMM_WORLD);

      /*Static Interleaved Scheduling to Send the data to corresponding process */
      if(myid = 0){
        for(i = 1;i < procs;i++){
          for(row = norm + 1 + i; row < N;row = row + procs){
            MPI_SEND(&A[row],N,MPI_float,i,0,MPI_COMM_WORLD);
            MPI_SEND(&B[row],1,MPI_float,i,0,MPI_COMM_WORLD);
          }
        }
         /*Elimination*/
        for (row = norm + 1; row < N; row = row + procs) {
          multiplier = A[row][norm] / A[norm][norm];
          for (col = norm; col < N; col++) {
            A[row][col] -= A[norm][col] * multiplier;
          }
          B[row] -= B[norm] * multiplier;
        }
        /*update data from processes*/
        for(i = 1;i < procs; i++){
          for (row = norm + 1 + i; row < N; row += procs) {
            MPI_Recv(&A[row],N,MPI_float,i,1,MPI_COMM_WORLD);
            MPI_Recv(&B[row],1,MPI_float,i,1,MPI_COMM_WORLD);
          }
        }
        if (norm == N - 2) {
          endTime = MPI_Wtime();
          printf("total elapsed time = %f\n", endTime - startTime);
        }
      }
      /* If not the processor 0, then receive data*/
    else{
      for (row = norm + 1 + myid; row < N; row += procs) {
        /*receive from process 0*/
        MPI_Recv(&A[row], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);   
        MPI_Recv(&B[row], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        /*Gaussian elimination*/
        multiplier = A[row][norm] / A[norm][norm];
        for (col = norm; col < N; col++) {
            A[row][col] -= A[norm][col] * multiplier;
          }
          B[row] -= B[norm] * multiplier;
          /*SEND BACK TO PROCESSOR 0*/
          MPI_SEND(&A[row],N,MPI_FLOAT,0,1,MPI_COMM_WORLD);
          MPI_SEND(&B[row],1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
    }

  /* (Diagonal elements are not normalized to 1.  This is treated in back
   * substitution.)
   */
  MPI_Barrier(MPI_COMM_WORLD);
  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];
    }
    X[row] /= A[row][row];
  }
}

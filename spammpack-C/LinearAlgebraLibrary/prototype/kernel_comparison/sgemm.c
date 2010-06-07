#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int
main (int argc, char **argv)
{
  int loops = 100;
  int N = 32;

  int i, j, loop;

  float alpha = 1.0;
  float beta = 1.0;
  float *A, *B, *C;

  struct timeval startBlas, stopBlas;
  double timeBlas;

  /* Read command line. */
  if (argc == 1)
  {
    /* No arguments. */
  }

  else if (argc == 3)
  {
    N = strtol(argv[1], NULL, 10);
    loops = strtol(argv[2], NULL, 10);
  }

  else
  {
    printf("wrong number of arguments\n");
    exit(1);
  }

  A = (float*) malloc(sizeof(float)*N*N);
  B = (float*) malloc(sizeof(float)*N*N);
  C = (float*) malloc(sizeof(float)*N*N);

  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
    {
      A[i+j*N] = rand()/(double) RAND_MAX;
      B[i+j*N] = rand()/(double) RAND_MAX;
      C[i+j*N] = 0;
    }
  }

  gettimeofday(&startBlas, NULL);
#ifdef DGEMM
  for (loop = 0; loop < loops; ++loop)
  {
    DGEMM("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);
  }
#else
#warning Need blas
#endif
  gettimeofday(&stopBlas, NULL);
  timeBlas = stopBlas.tv_sec-startBlas.tv_sec+(stopBlas.tv_usec-startBlas.tv_usec)/1.0e6;
  printf("%ix%i matrix (%i loop)\n", N, N, loops);
  printf("blas time: %f s total = %e s per loop iteration\n", timeBlas, timeBlas/(double) loops);

  free(A);
  free(B);
  free(C);
}

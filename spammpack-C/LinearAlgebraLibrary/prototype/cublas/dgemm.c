#include "config.h"

#ifdef HAVE_CUDA

#include <cublas.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//#define PRINTOUT__

int
main (int argc, char **argv)
{
  int loops = 100;
  int N = 100;
  float threshold = 1e-5;

  int i, j, loop;

  struct timeval startDMA1, stopDMA1;
  struct timeval startDMA2, stopDMA2;
  struct timeval startCublas, stopCublas;
  struct timeval startBlas, stopBlas;

  double timeDMA1, timeDMA2, timeCublas, timeBlas;

  float *A;
  float *B;
  float *C_cublas;
  float *C_blas;

  float alpha = 1.0;
  float beta = 0.0;

  float *d_A;
  float *d_B;
  float *d_C;

  float max_diff;
  int max_i, max_j;

  cublasStatus status;

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

  /* Allocate memory. */
  A = (float*) malloc(sizeof(float)*N*N);
  B = (float*) malloc(sizeof(float)*N*N);
  C_cublas = (float*) malloc(sizeof(float)*N*N);
  C_blas = (float*) malloc(sizeof(float)*N*N);

  /* Set A and B. */
  for (j = 0; j < N; ++j)
    for (i = 0; i < N; ++i) {
    {
      A[i+j*N] = rand()/(float) RAND_MAX;
      B[i+j*N] = rand()/(float) RAND_MAX;
      C_cublas[i+j*N] = 0.0;
      C_blas[i+j*N] = 0.0;
    }
  }

#ifdef PRINTOUT__
  printf("A =\n");
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" %f", A[i+j*N]);
    }
    printf("\n");
  }

  printf("B =\n");
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" %f", B[i+j*N]);
    }
    printf("\n");
  }
#endif

  printf("starting...\n");
  fflush(stdout);

  cublasInit();

  status = cublasAlloc(N*N, sizeof(float), (void**)&d_A);
  status = cublasAlloc(N*N, sizeof(float), (void**)&d_B);
  status = cublasAlloc(N*N, sizeof(float), (void**)&d_C);

  printf("timing DMA1...\n");
  fflush(stdout);
  gettimeofday(&startDMA1, NULL);
  for (loop = 0; loop < loops; ++loop)
  {
    cublasSetMatrix(N, N, sizeof(float), (void*) A, N, (void*) d_A, N);
    cublasSetMatrix(N, N, sizeof(float), (void*) B, N, (void*) d_B, N);
    cublasSetMatrix(N, N, sizeof(float), (void*) C_cublas, N, (void*) d_C, N);
    cudaThreadSynchronize();
  }
  gettimeofday(&stopDMA1, NULL);

  /* Multiply. */
  printf("timing cublas...\n");
  fflush(stdout);
  gettimeofday(&startCublas, NULL);
  for (loop = 0; loop < loops; ++loop)
  {
    cublasSgemm('N', 'N', N, N, N, alpha, d_A, N, d_B, N, beta, d_C, N);
    cudaThreadSynchronize();
  }
  gettimeofday(&stopCublas, NULL);

  printf("timing DMA2...\n");
  fflush(stdout);
  gettimeofday(&startDMA2, NULL);
  for (loop = 0; loop < loops; ++loop)
  {
    cublasGetMatrix(N, N, sizeof(float), (void*) d_C, N, (void*) C_cublas, N);
    cudaThreadSynchronize();
  }
  gettimeofday(&stopDMA2, NULL);

  cublasFree(d_A);
  cublasFree(d_B);
  cublasFree(d_C);
  cublasShutdown();

  printf("timing blas...\n");
  fflush(stdout);
  gettimeofday(&startBlas, NULL);
  for (loop = 0; loop < loops; ++loop)
  {
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C_blas, &N);
  }
  gettimeofday(&stopBlas, NULL);

  printf("comparing...\n");
  fflush(stdout);

  max_diff = 0.0;
  for (j = 0; j < N; ++j)
    for (i = 0; i < N; ++i) {
    {
      if (fabs(C_cublas[i+j*N]-C_blas[i+j*N]) > max_diff)
      {
        max_diff = fabs(C_cublas[i+j*N]-C_blas[i+j*N]);
        max_i = i;
        max_j = j;
      }
    }
  }

  if (max_diff > threshold)
  {
    printf("mismatch: C_cublas[%i][%i], cublas = %e, blas = %e, diff = %e\n",
        max_i, max_j, C_cublas[max_i+max_j*N], C_blas[max_i+max_j*N],
        fabs(C_cublas[max_i+max_j*N]-C_blas[max_i+max_j*N]));
  }

  timeDMA1 = stopDMA1.tv_sec-startDMA1.tv_sec+(stopDMA1.tv_usec-startDMA1.tv_usec)/1.0e6;
  timeDMA2 = stopDMA2.tv_sec-startDMA2.tv_sec+(stopDMA2.tv_usec-startDMA2.tv_usec)/1.0e6;
  timeCublas = stopCublas.tv_sec-startCublas.tv_sec+(stopCublas.tv_usec-startCublas.tv_usec)/1.0e6;
  timeBlas = stopBlas.tv_sec-startBlas.tv_sec+(stopBlas.tv_usec-startBlas.tv_usec)/1.0e6;

  printf("%ix%i matrix (%i loop)\n", N, N, loops);
  printf("DMA1 time: %f s total = %e s per loop iteration\n", timeDMA1, timeDMA1/(double) loops);
  printf("DMA2 time: %f s total = %e s per loop iteration\n", timeDMA2, timeDMA2/(double) loops);
  printf("cublas time: %f s total = %e s per loop iteration\n", timeCublas, timeCublas/(double) loops);
  printf("blas time: %f s total = %e s per loop iteration\n", timeBlas, timeBlas/(double) loops);
  printf("%i %e %e %e %e\n", N, timeDMA1/(double) loops, timeDMA2/(double) loops, timeCublas/(double) loops, timeBlas/(double) loops);

#ifdef PRINTOUT__
  printf("C_cublas =\n");
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" %f", C_cublas[i+j*N]);
    }
    printf("\n");
  }

  printf("C_blas =\n");
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" %f", C_blas[i+j*N]);
    }
    printf("\n");
  }
#endif
}

#else

int
main (int argc, char **argv)
{
}

#endif

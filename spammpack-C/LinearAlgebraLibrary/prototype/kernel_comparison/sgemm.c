#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>

int
main (int argc, char **argv)
{
  int loops = 100;
  unsigned int N = 32;

  int i, j, loop;

  float alpha = 1.2;
  float beta = 0.5;
  float *A, *B, *C;

  struct timeval start_blas, stop_blas;
  double walltime_blas, flops_blas;

  int parse;
  int longindex;
  char *short_options = "hN:l:";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
    { NULL, 0, NULL, 0 }
  };

  /* Read command line. */
  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h           This help\n");
        printf("-N N         Use NxN matrix blocks\n");
        printf("--loops N    Repeat each multiply N times\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

      default:
        printf("unknown command line argument\n");
        return -1;
        break;
    }
  }

  A = (float*) malloc(sizeof(float)*N*N);
  B = (float*) malloc(sizeof(float)*N*N);
  C = (float*) malloc(sizeof(float)*N*N);

  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
    {
      A[i+j*N] = rand()/(double) RAND_MAX;
      B[i+j*N] = rand()/(double) RAND_MAX;
      C[i+j*N] = rand()/(double) RAND_MAX;
    }
  }

  gettimeofday(&start_blas, NULL);
  for (loop = 0; loop < loops; ++loop)
  {
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);
  }
  gettimeofday(&stop_blas, NULL);
  walltime_blas = (stop_blas.tv_sec-start_blas.tv_sec+(stop_blas.tv_usec-start_blas.tv_usec)/1.0e6)/loops;
  flops_blas = ((double) N)*((double) N)*(2.*N+1.)/walltime_blas;
  printf("%ix%i matrix (%i loop)\n", N, N, loops);
  if (flops_blas < 1000*1000*1000)
  {
    printf("performance: total walltime = %f s, walltime/iteration = %e s = %1.2f Mflop/s\n", walltime_blas*loops, walltime_blas, flops_blas/1000./1000.);
  }

  else
  {
    printf("performance: total walltime = %f s, walltime/iteration = %e s = %1.2f Gflop/s\n", walltime_blas*loops, walltime_blas, flops_blas/1000./1000./1000.);
  }

  free(A);
  free(B);
  free(C);
}

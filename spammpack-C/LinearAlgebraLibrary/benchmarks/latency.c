#include <config.h>

#ifdef HAVE_CUDA

#include <cublas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

int
main (int argc, char **argv)
{
  float *A;
  void *device_pointer;
  struct timeval start, stop;
  double time_elapsed;

  int i, j;
  int N = 250;
  int number_tests = 1;

  if (argc == 3)
  {
    N = strtol(argv[1], NULL, 10);
    number_tests = strtol(argv[2], NULL, 10);
  }

  A = (float*) malloc(sizeof(float)*N*N);
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      A[i+j*N] = rand()/(float) RAND_MAX;
    }
  }

  cublasAlloc(N*N, sizeof(float), &device_pointer);

  printf("%ix%i matrix, running %i tests\n", N, N, number_tests);
  gettimeofday(&start, NULL);
  for (i = 0; i < number_tests; ++i)
  {
    cublasSetMatrix(N, N, sizeof(float), (void*) A, N, device_pointer, N);
  }
  gettimeofday(&stop, NULL);

  cublasFree(device_pointer);

  time_elapsed = stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6;

  printf("time elapsed: %e s\n", time_elapsed);
  printf("time elapsed per iteration: %e s\n", time_elapsed/(double) number_tests);

  return 0;
}

#else

int
main ()
{
  return 0;
}

#endif

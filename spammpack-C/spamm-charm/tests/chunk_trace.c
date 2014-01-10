#include "chunk.h"

#include <math.h>
#include <stdio.h>

int
main (int argc, char **argv)
{
  const int N = 512;
  const int N_chunk = 512;
  const int N_basic = 4;

  void *A = chunk_alloc(N_chunk, N_basic, N, 0, 0);

  double *A_dense = calloc(N_chunk*N_chunk, sizeof(double));

  for(int i = 0; i < N_chunk*N_chunk; i++)
  {
    A_dense[i] = rand()/(double) RAND_MAX;
  }

  chunk_set(A, A_dense);

  double trace = chunk_trace(A);

  printf("done calculating trace, verifying...\n");

  double trace_exact = 0;
  for(int i = 0; i < N_chunk; i++)
  {
    trace_exact += A_dense[i*N_chunk+i];
  }

  if(fabs(trace_exact-trace) > 1e-10)
  {
    printf("mismatch\n");
    return -1;
  }

  printf("matrices are identical\n");

  free(A_dense);
  free(A);

  return 0;
}

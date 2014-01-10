#include "config.h"

#include "chunk.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

int
main (int argc, char **argv)
{
  const int N = 2048;
  const int N_chunk = 2048;
  const int N_basic = 4;

  void *A = chunk_alloc(N_chunk, N_basic, N, 0, 0);
  void *C = chunk_alloc(N_chunk, N_basic, N, 0, 0);

  double *A_dense = calloc(N_chunk*N_chunk, sizeof(double));

  for(int i = 0; i < N_chunk*N_chunk; i++)
  {
    A_dense[i] = rand()/(double) RAND_MAX;
  }

  chunk_set(A, A_dense);

  struct timespec start_time;
  clock_gettime(CLOCKTYPE, &start_time);
  chunk_add(0.0, C, 1.2, A);
  struct timespec end_time;
  clock_gettime(CLOCKTYPE, &end_time);

  printf("done adding chunks, %1.2f seconds, verifying...\n",
      (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
      (start_time.tv_sec+start_time.tv_nsec/1.0e9));

  double *C_dense = chunk_to_dense(C);

  for(int i = 0; i < N_chunk; i++)
  {
    for(int j = 0; j < N_chunk; j++)
    {
      double C_exact = 1.2*A_dense[i*N_chunk+j];

      if(fabs(C_exact-C_dense[i*N_chunk+j]) > 1e-10)
      {
        printf("mismatch C[%d][%d]\n", i, j);
        return -1;
      }
    }
  }

  printf("matrices are identical\n");

  free(A_dense);
  free(C_dense);
  free(A);
  free(C);

  return 0;
}

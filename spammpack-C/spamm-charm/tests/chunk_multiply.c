#include "config.h"

#include "chunk.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

int
main (int argc, char **argv)
{
  const int N = 900;
  const int N_chunk = 900;
  const int N_basic = 5;

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
  chunk_multiply(0.0, A, A, C);
  struct timespec end_time;
  clock_gettime(CLOCKTYPE, &end_time);

  printf("done multiplying chunks, %1.2f seconds, verifying...\n",
      (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
      (start_time.tv_sec+start_time.tv_nsec/1.0e9));

  double *C_dense = chunk_to_dense(C);

  double C_exact;
  for(int i = 0; i < N_chunk; i++)
  {
    for(int j = 0; j < N_chunk; j++)
    {
      C_exact = 0;
      for(int k = 0; k < N_chunk; k++)
      {
        C_exact += A_dense[i*N_chunk+k]*A_dense[k*N_chunk+j];
      }

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

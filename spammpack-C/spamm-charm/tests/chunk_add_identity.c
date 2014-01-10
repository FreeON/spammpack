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

  chunk_add_identity(1.2, A);

  printf("done adding to matrix, verifying...\n");

  double *C_dense = chunk_to_dense(A);

  for(int i = 0; i < N_chunk; i++)
  {
    for(int j = 0; j < N_chunk; j++)
    {
      double C_exact = A_dense[i*N_chunk+j];
      if(i == j)
      {
        C_exact += 1.2;
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

  return 0;
}

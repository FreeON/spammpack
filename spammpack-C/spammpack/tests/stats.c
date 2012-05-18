#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

#define N 1000
#define FILL 0.60

int
main ()
{
  int result = 0;

  unsigned int i, j;
  unsigned int nonzeros = 0;

  struct spamm_hashed_t *A;
  float *A_dense;

  A_dense = (float*) malloc(sizeof(float)*N*N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (rand()/(double) RAND_MAX > FILL)
      {
        A_dense[i*N+j] = 1.0;
        nonzeros++;
      }

      else
      {
        A_dense[i*N+j] = 0.0;
      }
    }
  }

  A = spamm_convert_dense_to_spamm(N, N, row_major, A_dense, row_major);
  if (spamm_number_nonzero(A) != nonzeros)
  {
    printf("found %u nonzeros, should have found %u\n", spamm_number_nonzero(A), nonzeros);
    result = -1;
  }

  spamm_hashed_delete(&A);
  free(A_dense);

  return result;
}

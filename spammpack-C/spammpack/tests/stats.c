#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

#define FILL 0.60

int
main ()
{
  int result = 0;

  unsigned int i, j;
  unsigned int nonzeros = 0;

  const unsigned int N[] = { 1000, 1000 };
  const unsigned int contiguous_tier = 5;
  const unsigned int N_block = 4;
  const short use_linear_tree = 1;

  struct spamm_matrix_t *A;
  float *A_dense;

  A_dense = (float*) malloc(sizeof(float)*N[0]*N[1]);
  for (i = 0; i < N[0]; i++) {
    for (j = 0; j < N[1]; j++)
    {
      if (rand()/(double) RAND_MAX > FILL)
      {
        A_dense[i*N[1]+j] = 1.0;
        nonzeros++;
      }

      else
      {
        A_dense[i*N[1]+j] = 0.0;
      }
    }
  }

  A = spamm_convert_dense_to_spamm(2, N, contiguous_tier, N_block, use_linear_tree, row_major, A_dense, row_major);
  if (spamm_number_nonzero(A) != nonzeros)
  {
    printf("found %u nonzeros, should have found %u\n", spamm_number_nonzero(A), nonzeros);
    result = -1;
  }

  spamm_delete(&A);
  free(A_dense);

  return result;
}

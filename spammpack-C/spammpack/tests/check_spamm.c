#include "spamm.h"
#include <stdlib.h>

int
main ()
{
  int result = 0;

  unsigned int N = 129;
  const unsigned int linear_tier = 2;
  const unsigned int contiguous_tier = 3;
  unsigned int i, j;

  struct spamm_matrix_t *A;

  A = spamm_new(N, N, linear_tier, contiguous_tier, row_major);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      spamm_set(i, j, rand()/(float) RAND_MAX, A);
    }
  }

  result = spamm_check(A, 1e-7);

  spamm_delete(&A);

  return result;
}

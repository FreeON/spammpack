#include "spamm.h"

int
main ()
{
  int result = 0;

  const unsigned int N = 1000;
  const unsigned int linear_tier = 4;
  const unsigned int contiguous_tier = 5;

  struct spamm_matrix_t *A;

  A = spamm_new(N, N, linear_tier, contiguous_tier, row_major);

  spamm_expand(A);

  spamm_delete(&A);

  return result;
}

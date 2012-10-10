#include "spamm.h"
#include <stdlib.h>

int
main ()
{
  int result = 0;

  unsigned int N[] = { 129, 129 };
  const unsigned int linear_tier = 2;
  const unsigned int contiguous_tier = 3;
  unsigned int i[2];

  struct spamm_matrix_t *A;

  A = spamm_new(2, N, linear_tier, contiguous_tier, row_major);

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      spamm_set(i, rand()/(float) RAND_MAX, A);
    }
  }

  result = spamm_check(A, 1e-7);

  spamm_delete(&A);

  return result;
}

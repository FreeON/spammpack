#include "spamm.h"
#include <stdlib.h>

int
main ()
{
  int result = 0;

  unsigned int N[] = { 129, 129 };
  const unsigned int contiguous_tier = 3;
  const short use_linear_tree = 1;
  unsigned int i[2];

  struct spamm_matrix_t *A;

  A = spamm_new(2, N, contiguous_tier, use_linear_tree);

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

#include "spamm.h"
#include <stdlib.h>

int
main ()
{
  int result = 0;

  unsigned int N = 129;
  unsigned int i, j;

  struct spamm_hashed_t *A;

  //A = spamm_hashed_new(N, N, row_major);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      spamm_hashed_set(i, j, rand()/(float) RAND_MAX, A);
    }
  }

  result = spamm_hashed_check(A, 1e-7);

  spamm_hashed_delete(&A);

  return result;
}

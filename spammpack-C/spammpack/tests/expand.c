#include "spamm.h"

#define N 1000

int
main ()
{
  int result = 0;
  struct spamm_hashed_t *A;

  A = spamm_hashed_new(N, N, row_major);

  spamm_expand(A);

  spamm_hashed_delete(&A);

  return result;
}

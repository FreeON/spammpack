#include "spamm.h"

#define N 1000

int
main ()
{
  int result = 0;
  struct spamm_t *A;

  A = spamm_new(N, N);

  spamm_expand(A);

  spamm_delete(&A);

  return result;
}

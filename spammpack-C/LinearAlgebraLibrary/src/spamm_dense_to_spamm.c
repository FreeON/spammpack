#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_dense_to_spamm (const double *A_dense, struct spamm_t *A)
{
  assert(A_dense != NULL);
  assert(A != NULL);

  int i, j;

  /* Allocate new tree. */
}

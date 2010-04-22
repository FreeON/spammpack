#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_print_spamm (const struct spamm_t *A)
{
  int i, j;

  assert(A != NULL);

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      printf(" % f", spamm_get(i, j, A));
    }
    printf("\n");
  }
}

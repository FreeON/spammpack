#include "spamm.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A)
{
  float Aij = 0;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    printf("illegal index values for A_ij\n");
    exit(1);
  }

  return Aij;
}

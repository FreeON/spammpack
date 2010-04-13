#include "spamm.h"
#include <stdlib.h>

void
spamm_new (struct spamm_t *A)
{
  A->M = 0;
  A->N = 0;
  A->root = NULL;
}

#include "spamm.h"
#include <stdlib.h>

void
spamm_delete (struct spamm_t **A)
{
  /* Delete all data on all tiers. */

  /* Delete the matrix. */
  free(*A);
  *A = NULL;
}

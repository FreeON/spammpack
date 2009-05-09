#include "lal.h"

void
lal_free (lal_matrix_t **A)
{
  free((*A)->data);
  free(*A);
}

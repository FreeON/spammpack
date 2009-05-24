#include "lal.h"

#include <stdlib.h>

void
f90_lal_set_ (int *i, int *j, double *Aij, int *f90_A)
{
  struct lal_matrix_t *A = NULL;

  lal_integer_to_pointer(f90_A, &A);

  lal_set(*i, *j, *Aij, A);
}

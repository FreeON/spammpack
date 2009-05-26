#include "lal.h"

#include <stdlib.h>

void
f90_lal_free_ (int *f90_A)
{
  struct lal_matrix_t *A = NULL;

  lal_integer_to_pointer(f90_A, &A);

  lal_free(&A);

  lal_pointer_to_integer(f90_A, A);
}

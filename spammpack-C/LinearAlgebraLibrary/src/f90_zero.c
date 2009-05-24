#include "lal.h"

#include <stdio.h>

void
f90_lal_zero_ (int *f90_A)
{
  struct lal_matrix_t *A = NULL;

  lal_integer_to_pointer(f90_A, &A);

  lal_zero(A);
}

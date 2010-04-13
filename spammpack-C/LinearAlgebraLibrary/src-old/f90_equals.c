#include "lal.h"

#include <stdlib.h>

int
f90_lal_equals_ (int *f90_A, int *f90_B, double *tolerance)
{
  struct lal_matrix_t *A = NULL;
  struct lal_matrix_t *B = NULL;

  lal_integer_to_pointer(f90_A, &A);
  lal_integer_to_pointer(f90_B, &B);

  return lal_equals(A, B, *tolerance);
}

#include "lal.h"

#include <stdlib.h>

int
f90_lal_allocate_ (int *M, int *N, int *f90_A)
{
  int return_value;

  struct lal_matrix_t *A = NULL;

  return_value = lal_allocate(*M, *N, &A);

  lal_pointer_to_integer(f90_A, A);

  return return_value;
}

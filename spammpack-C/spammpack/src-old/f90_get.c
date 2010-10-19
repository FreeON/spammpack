#include "lal.h"

#include <stdlib.h>

double
f90_lal_get_ (int *i, int *j, int * f90_A)
{
  struct lal_matrix_t *A = NULL;

  lal_integer_to_pointer(f90_A, &A);

  return lal_get((*i)-1, (*j)-1, A);
}

#include "lal.h"

void
lal_pointer_to_integer (int *f90_A, struct lal_matrix_t *A)
{
  int i;
  char *pointer;

  pointer = (char*) (&A);

  for (i = 0; i < sizeof(void *); ++i)
  {
    f90_A[i] = (int) pointer[i];
  }
}

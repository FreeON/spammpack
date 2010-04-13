#include "lal.h"

void
lal_integer_to_pointer (int *f90_A, struct lal_matrix_t **A)
{
  int i;
  char *pointer;
  char temp;

  pointer = (char*) A;

  for (i = 0; i < sizeof(void *); ++i)
  {
    pointer[i] = (char) f90_A[i];
  }
}

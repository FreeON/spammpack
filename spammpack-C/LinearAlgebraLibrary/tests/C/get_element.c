#include <lal.h>

#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  lal_matrix_t *A = NULL;

  lal_allocate(10, 10, &A);

  lal_zero(A);

  if (lal_get(0, 0, A) == 0.0) { return 0; }
  else { return 1; }
}

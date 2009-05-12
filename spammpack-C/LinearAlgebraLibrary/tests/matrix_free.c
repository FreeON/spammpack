#include <lal.h>

int
main ()
{
  lal_matrix_t *A;

  lal_allocate(100, 100, &A);

  lal_free(&A);

  return 0;
}

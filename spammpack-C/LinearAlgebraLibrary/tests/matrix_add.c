#include <lal.h>

int
main ()
{
  int M = 10;
  int N = 10;

  lal_matrix_t *A;
  lal_matrix_t *B;
  lal_matrix_t *C;

  lal_allocate(M, N, &A);
  lal_allocate(M, N, &B);
  lal_allocate(M, N, &C);

  lal_rand(A);
  lal_rand(B);

  lal_add(A, B, C);

  return 0;
}

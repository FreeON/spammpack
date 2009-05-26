#include <lal.h>

#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int M = 4;
  int N = 4;

  int i, j, k;

  lal_matrix_t *A = NULL;
  lal_matrix_t *B = NULL;
  lal_matrix_t *C = NULL;
  lal_matrix_t *C_reference = NULL;

  if (lal_allocate(M, N, &A) != 0) { return 1; }
  if (lal_allocate(M, N, &B) != 0) { return 1; }
  if (lal_allocate(M, N, &C) != 0) { return 1; }
  if (lal_allocate(M, N, &C_reference) != 0) { return 1; }

  lal_rand(A);
  lal_rand(B);

  /* Multiply by hand. */
  lal_zero(C_reference);
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < M; ++k)
      {
        lal_set(i, j, lal_get(i, j, C_reference)+lal_get(i, k, A)*lal_get(k, j, B), C_reference);
      }
    }
  }

  /* Multiply by library. */
  lal_dgemm("N", "N", M, N, N, 1.0, A, N, B, M, 1.0, C, N);

  if (lal_equals(C, C_reference) != 0)
  {
    printf("[matrix_multiply] C is not equal to C_reference\n");
    printf("[matrix_multiply] C =\n");
    lal_print(C);
    printf("[matrix_multiply] C_reference =\n");
    lal_print(C_reference);
    return -1;
  }

  else { return 0; }
}

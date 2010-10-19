#include "lal.h"

#include <stdlib.h>

void
f90_lal_dgemm_ (char *transA, char *transB, int *M, int *N,
    int *K, double *alpha, int *f90_A, int *lda,
    int *f90_B, int *ldb, double *beta, int *f90_C,
    int *ldc)
{
  struct lal_matrix_t *A;
  struct lal_matrix_t *B;
  struct lal_matrix_t *C;

  lal_integer_to_pointer(f90_A, &A);
  lal_integer_to_pointer(f90_B, &B);
  lal_integer_to_pointer(f90_C, &C);

  lal_dgemm(transA, transB, *M, *N, *K, *alpha, A, *lda, B, *ldb, *beta, C, *ldc);
}

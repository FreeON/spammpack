#include "lal.h"

void
lal_dgemm_ (const char *transA, const char *transB, const int M, const int N,
    const int K, const double alpha, const lal_matrix_t *A, const int lda,
    const lal_matrix_t *B, const int ldb, const double beta, lal_matrix_t *C,
    const int ldc)
{
  lal_dgemm(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

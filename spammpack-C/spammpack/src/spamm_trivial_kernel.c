#include "spamm.h"

/** Trivial sgemm kernel.
 *
 * This kernel multiplies nodes of a matrix tree.
 *
 * \f$ C = \alpha A \times B \f$
 *
 * \bug Does not support all operations yet. Currently only opA == 'N' and opB
 * == 'N' are supported.
 *
 * @param[in] opA opA specifies the form of op( A ) to be used in the matrix
 * multiplication as follows:
 * - OpA = 'N' or 'n', op(A) = A.
 * - OpA = 'T' or 't', op(A) = A'.
 * - OpA = 'C' or 'c', op(A) = A'.
 * @param[in] opB opA specifies the form of op( B ) to be used in the matrix
 * multiplication as follows:
 * - opB = 'N' or 'n', op(B) = B.
 * - opB = 'T' or 't', op(B) = B'.
 * - opB = 'C' or 'c', op(B) = B'.
 * @param[in] M M specifies the number of rows of the matrix op(A) and of
 * the matrix C. M must be at least zero.
 * @param[in] N N specifies the number of columns of the matrix op(B) and the
 * number of columns of the matrix C. N must be at least zero.
 * @param[in] K K specifies the number of columns of the matrix op(A) and the
 * number of rows of the matrix op(B). K must be at least zero.
 * @param[in] alpha alpha specifies the scalar alpha.
 * @param[in] A_block_dense Array of DIMENSION (lda, Ka), where ka is K when
 * opA = 'N' or 'n', and is M otherwise. Before entry with opA = 'N' or 'n',
 * the leading M by K part of the array A must contain the matrix A, otherwise
 * the leading K by M part of the array A must contain the matrix A.
 * @param[in] lda lda specifies the first dimension of A as declared in the
 * calling (sub) program. When opA = 'N' or 'n' then lda must be at least
 * max(1, M), otherwise lda must be at least max(1, K).
 * @param[in] B_block_dense Array of DIMENSION (ldb, Kb), where Kb is N when
 * opB = 'N' or 'n', and is K otherwise. Before entry with opB = 'N' or 'n',
 * the leading K by N part of the array B must contain the matrix B, otherwise
 * the leading N by K part of the array B must contain the matrix B.
 * @param[in] ldb ldb specifies the first dimension of B as declared in the
 * calling (sub) program. When opB = 'N' or 'n' then ldb must be at least
 * max(1, K), otherwise ldb must be at least max(1, N).
 * @param[in] beta beta specifies the scalar beta. When beta is supplied
 * as zero then C need not be set on input.
 * @param[out] C_block_dense Array of DIMENSION (ldc, N). Before entry, the
 * leading m by n part of the array C must contain the matrix C, except when
 * beta is zero, in which case C need not be set on entry. On exit, the array
 * C is overwritten by the M by N matrix (alpha*op(A)*op(B) + beta*C).
 * @param[in] ldc ldc specifies the first dimension of C as declared in the
 * calling (sub) program. ldc must be at least max(1, M).
 */
void
spamm_sgemm_trivial (const char opA, const char opB,
    const unsigned int M, const unsigned int N, const unsigned int K,
    const floating_point_t alpha,
    const floating_point_t *A_block_dense,
    const unsigned int lda,
    const floating_point_t *B_block_dense,
    const unsigned int ldb,
    const floating_point_t beta,
    floating_point_t *C_block_dense,
    const unsigned int ldc)
{
  unsigned int i, j, k;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < K; k++)
      {
        C_block_dense[spamm_dense_index(i, j, M, N)] += alpha*A_block_dense[spamm_dense_index(i, k, M, K)]*B_block_dense[spamm_dense_index(k, j, K, N)];
      }
    }
  }
}

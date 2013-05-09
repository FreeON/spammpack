/** @file */

#include "spamm.h"

#define CONCAT2(a, b) a ## _ ## b
#define FUNC(a, b) CONCAT2(a, b)

/** A recursive implementation of {s,d}gemm(). This function is not feature
 * complete, it hardly does anything, but multiply two matrices.
 */
void FUNC(spamm, FUNCNAME) (char * transA, char * transB,
    int *M, int *N, int *K,
    FUNCTYPE *alpha, FUNCTYPE *A, int *LDA, FUNCTYPE *B, int *LDB,
    FUNCTYPE *beta, FUNCTYPE *C, int *LDC)
{
  int i, j, k;

  if(*transA != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  if(*transB != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  for(i = 0; i < *M; i++) {
    for(j = 0; j < *N; j++)
    {
      C[spamm_index_column_major(i, j, *M, *N)] *= (*beta);
      for(k = 0; k < *K; k++)
      {
        C[spamm_index_column_major(i, j, *M, *N)] += (*alpha)*A[spamm_index_column_major(i, k, *M, *K)]*B[spamm_index_column_major(k, j, *K, *N)];
      }
    }
  }
}

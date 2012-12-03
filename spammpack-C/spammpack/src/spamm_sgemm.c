/** @file */

#include "spamm.h"

/** A recursive implementation of sgemm(). This function is not feature complete,
 * it hardly does anything.
 */
void spamm_sgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    float *alpha, float *A, int *LDA, float *B, int *LDB,
    float *beta, float *C, int *LDC)
{
  int i, j, k;

  SPAMM_WARN("using SpAMM sgemm()\n");

  if(*transA != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  if(*transB != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  for(i = 0; i < *M; i++) {
    for(j = 0; j < *N; j++) {
      for(k = 0; k < *K; k++)
      {
        C[spamm_index_column_major(i, j, *M, *N)] = (*beta)*C[spamm_index_column_major(i, j, *M, *N)]
          +(*alpha)*A[spamm_index_column_major(i, k, *M, *K)]*B[spamm_index_column_major(k, j, *K, *N)];
      }
    }
  }
}

/** A recursive implementation of dgemm(). This function is not feature complete,
 * it hardly does anything.
 */
void spamm_dgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    double *alpha, double *A, int *LDA, double *B, int *LDB,
    double *beta, double *C, int *LDC)
{
  int i, j, k;

  SPAMM_WARN("using SpAMM dgemm()\n");

  if(*transA != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  if(*transB != 'N')
  {
    SPAMM_FATAL("FIXME\n");
  }

  for(i = 0; i < *M; i++) {
    for(j = 0; j < *N; j++) {
      for(k = 0; k < *K; k++)
      {
        C[spamm_index_column_major(i, j, *M, *N)] = (*beta)*C[spamm_index_column_major(i, j, *M, *N)]
          +(*alpha)*A[spamm_index_column_major(i, k, *M, *K)]*B[spamm_index_column_major(k, j, *K, *N)];
      }
    }
  }
}

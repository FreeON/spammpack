/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "chunk_block.h"

#include <assert.h>

/** Calculate a row-major offset. */
#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))

/** A basic matrix product. The matrix elements are assumed to be in row-major
 * order.
 *
 * @param A The pointer to matrix A.
 * @param B The pointer to matrix B.
 * @param C The pointer to matrix C.
 * @param N The size of the matrices.
 */
void
chunk_block_multiply (const double *const restrict A,
    const double *const restrict B,
    double *const restrict C,
    const int N)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

#ifdef CHUNK_BLOCK_TRANSPOSE
  if(N > CHUNK_BLOCK_TRANSPOSE)
  {
    double *const restrict B_transpose = calloc(N*N, sizeof(double));
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        B_transpose[ROW_MAJOR(i, j, N)] = B[ROW_MAJOR(j, i, N)];
      }
    }

    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          C[ROW_MAJOR(i, j, N)] += A[ROW_MAJOR(i, k, N)] * B_transpose[ROW_MAJOR(j, k, N)];
        }
      }
    }

    free(B_transpose);
  }

  else
  {
#endif
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          C[ROW_MAJOR(i, j, N)] += A[ROW_MAJOR(i, k, N)] * B[ROW_MAJOR(k, j, N)];
        }
      }
    }
#ifdef CHUNK_BLOCK_TRANSPOSE
  }
#endif
}

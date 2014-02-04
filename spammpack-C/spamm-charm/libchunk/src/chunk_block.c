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
#if defined(BLOCK_SLEEP) || defined(CHUNK_BLOCK_NO_WORK)
#include <time.h>
#endif

/** Calculate a row-major offset. */
#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))

/** Calculate a column-major offset. */
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))

/** A basic matrix product. The matrix elements are assumed to be in row-major
 * order.
 *
 * @param A The pointer to matrix A.
 * @param B The pointer to matrix B.
 * @param C The pointer to matrix C.
 * @param N The size of the matrices.
 */
void
chunk_block_multiply_general (const double *const restrict A,
    const double *const restrict B,
    double *const restrict C,
    const int N)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

#ifdef CHUNK_BLOCK_NO_WORK
  struct timespec requested_sleep;
  requested_sleep.tv_sec = 0;
  requested_sleep.tv_nsec = 5e6;

  nanosleep(&requested_sleep, NULL);
#else

#ifdef CHUNK_BLOCK_TRANSPOSE
  if(N >= CHUNK_BLOCK_TRANSPOSE)
  {
    double *const restrict B_transpose = calloc(N*N, sizeof(double));
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        B_transpose[COLUMN_MAJOR(i, j, N)] = B[COLUMN_MAJOR(j, i, N)];
      }
    }

    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          C[COLUMN_MAJOR(i, j, N)] += A[COLUMN_MAJOR(i, k, N)] * B_transpose[COLUMN_MAJOR(j, k, N)];
        }
      }
    }

    free(B_transpose);

#ifdef BLOCK_SLEEP
    /* Sleep to add some extra slowness. */
    struct timespec requested_sleep;
    requested_sleep.tv_sec = 0;
    requested_sleep.tv_nsec = BLOCK_SLEEP;

    nanosleep(&requested_sleep, NULL);
#endif
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
          C[COLUMN_MAJOR(i, j, N)] += A[COLUMN_MAJOR(i, k, N)] * B[COLUMN_MAJOR(k, j, N)];
        }
      }
    }
#ifdef CHUNK_BLOCK_TRANSPOSE
  }
#endif

#endif
}

/** A basic matrix product. The matrix elements are assumed to be in row-major
 * order.
 *
 * @param A The pointer to matrix A.
 * @param B The pointer to matrix B.
 * @param C The pointer to matrix C.
 */
void
chunk_block_multiply_4_x_4 (const double *const restrict A,
    const double *const restrict B,
    double *const restrict C)
{
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume_aligned(C, 64);

  double A_dilated[4*4*4] __attribute__((aligned(64)));

  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      A_dilated[COLUMN_MAJOR(i, j, 4)+0] = A[COLUMN_MAJOR(i, j, 4)];
      A_dilated[COLUMN_MAJOR(i, j, 4)+1] = A[COLUMN_MAJOR(i, j, 4)];
      A_dilated[COLUMN_MAJOR(i, j, 4)+2] = A[COLUMN_MAJOR(i, j, 4)];
      A_dilated[COLUMN_MAJOR(i, j, 4)+3] = A[COLUMN_MAJOR(i, j, 4)];
    }
  }

  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        C[COLUMN_MAJOR(i, j, 4)] += A[COLUMN_MAJOR(i, k, 4)] * B[COLUMN_MAJOR(k, j, 4)];
      }
    }
  }
}

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
  //chunk_block_multiply_general(A, B, C, N);
  chunk_block_multiply_4_x_4(A, B, C);
}

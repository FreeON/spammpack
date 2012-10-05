#include "spamm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_TOLERANCE 5e-7
#define DECAY 0.7

//#define PRINT_MATRIX

int
main (int argc, char **argv)
{
  int result;

  const unsigned int N = 513;
  const unsigned int linear_tier = 3;
  const unsigned int contiguous_tier = 4;

  double alpha = 1.2;
  double beta = 0.8;

  float *A_dense;
  float *B_dense;

  float max_diff;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;

  int i, j;

  short sparse_test;
  const short sparse_A[] = { 0, 0, 1, 1 };
  const short sparse_B[] = { 0, 1, 0, 1 };

  for (sparse_test = 0; sparse_test < 4; sparse_test++)
  {
    A_dense = calloc(N*N, sizeof(double));
    B_dense = calloc(N*N, sizeof(double));

    for (i = 0; i < N; i++)
    {
      A_dense[spamm_index_row_major(i, i, N, N)] = 0.8+rand()/(double) RAND_MAX;
      B_dense[spamm_index_row_major(i, i, N, N)] = 0.8+rand()/(double) RAND_MAX;

      for (j = i+1; j < N; j++)
      {
        /* Exponential decay. */
        if (sparse_A[sparse_test] == 1)
        {
          A_dense[spamm_index_row_major(i, j, N, N)] = A_dense[spamm_index_row_major(i, i, N, N)]*exp(-fabsf(i-j)/DECAY);
          A_dense[spamm_index_row_major(j, i, N, N)] = A_dense[spamm_index_row_major(i, j, N, N)];
        }

        else
        {
          A_dense[spamm_index_row_major(i, j, N, N)] = 0.8+rand()/(double) RAND_MAX;
          A_dense[spamm_index_row_major(j, i, N, N)] = 0.8+rand()/(double) RAND_MAX;
        }

        if (sparse_B[sparse_test] == 1)
        {
          B_dense[spamm_index_row_major(i, j, N, N)] = B_dense[spamm_index_row_major(i, i, N, N)]*exp(-fabsf(i-j)/DECAY);
          B_dense[spamm_index_row_major(j, i, N, N)] = B_dense[spamm_index_row_major(i, j, N, N)];
        }

        else
        {
          B_dense[spamm_index_row_major(i, j, N, N)] = 0.8+rand()/(double) RAND_MAX;
          B_dense[spamm_index_row_major(j, i, N, N)] = 0.8+rand()/(double) RAND_MAX;
        }
      }
    }

#ifdef PRINT_MATRIX
    printf("A_dense =\n");
    spamm_print_dense(N, N, row_major, A_dense);

    printf("B_dense =\n");
    spamm_print_dense(N, N, row_major, B_dense);
#endif

    /* Convert to SpAMM matrix. */
    A = spamm_convert_dense_to_spamm(N, N, linear_tier, contiguous_tier, row_major, A_dense, row_major);
    B = spamm_convert_dense_to_spamm(N, N, linear_tier, contiguous_tier, row_major, B_dense, row_major);

    /* Add by hand for verification. */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
      {
        A_dense[spamm_index_row_major(i, j, N, N)] =
          alpha*A_dense[spamm_index_row_major(i, j, N, N)]
          + beta*B_dense[spamm_index_row_major(i, j, N, N)];
      }
    }

#ifdef PRINT_MATRIX
    printf("A_dense =\n");
    spamm_print_dense(N, N, row_major, A_dense);
#endif

    spamm_add(alpha, A, beta, B);

#ifdef PRINT_MATRIX
    printf("A =\n");
    spamm_print(A);
#endif

    /* Check tree consistency. */
    result = spamm_check(A, TEST_TOLERANCE);

    /* Compare result. */
    max_diff = 0.0;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
      {
        if (fabs(A_dense[spamm_index_row_major(i, j, N, N)]-spamm_get(i, j, A)) > max_diff)
        {
          max_diff = fabs(A_dense[spamm_index_row_major(i, j, N, N)]-spamm_get(i, j, A));
        }
      }
    }

    /* Test result. */
    if (max_diff > TEST_TOLERANCE)
    {
      result = SPAMM_ERROR;
      printf("[add_spamm] max diff (sparse_A = %u, sparse_B = %u, test tolerance was %1.2e) = %e\n",
          sparse_A[sparse_test], sparse_B[sparse_test], TEST_TOLERANCE, max_diff);
    }

    free(A_dense);
    free(B_dense);

    spamm_delete(&A);
    spamm_delete(&B);
  }

  return result;
}

#include "spamm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_TOLERANCE 5e-7
#define DECAY 2.5

//#define PRINT_MATRIX

int
main (int argc, char **argv)
{
  int result;

  const unsigned int N[] = { 734, 734 };
  const unsigned int chunk_tier = 4;
  const short use_linear_tree = 0;

  double alpha = 1.2;
  double beta = 0.8;

  float *A_dense;
  float *B_dense;

  float max_diff;
  unsigned int max_diff_i[2];

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;

  unsigned int i[2];

  short sparse_test;
  const short sparse_A[] = { 0, 0, 1, 1 };
  const short sparse_B[] = { 0, 1, 0, 1 };

  double flop;

  for(sparse_test = 0; sparse_test < 4; sparse_test++)
  {
    printf("sparse_A = %i, sparse_B = %i, ", sparse_A[sparse_test], sparse_B[sparse_test]);

    A_dense = calloc(N[0]*N[1], sizeof(double));
    B_dense = calloc(N[0]*N[1], sizeof(double));

    for(i[0] = 0; i[0] < N[0]; i[0]++)
    {
      A_dense[spamm_index_row_major(i[0], i[0], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;
      B_dense[spamm_index_row_major(i[0], i[0], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;

      for(i[1] = i[0]+1; i[1] < N[1]; i[1]++)
      {
        /* Exponential decay. */
        if(sparse_A[sparse_test] == 1)
        {
          A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = A_dense[spamm_index_row_major(i[0], i[0], N[0], N[1])]*exp(-fabsf(i[0]-i[1])/DECAY);
          A_dense[spamm_index_row_major(i[1], i[0], N[0], N[1])] = A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])];
        }

        else
        {
          A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;
          A_dense[spamm_index_row_major(i[1], i[0], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;
        }

        if(sparse_B[sparse_test] == 1)
        {
          B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = B_dense[spamm_index_row_major(i[0], i[0], N[0], N[1])]*exp(-fabsf(i[0]-i[1])/DECAY);
          B_dense[spamm_index_row_major(i[1], i[0], N[0], N[1])] = B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])];
        }

        else
        {
          B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;
          B_dense[spamm_index_row_major(i[1], i[0], N[0], N[1])] = 0.8+rand()/(double) RAND_MAX;
        }
      }
    }

#ifdef PRINT_MATRIX
    printf("A_dense =\n");
    spamm_print_dense(N[0], N[1], row_major, A_dense);

    printf("B_dense =\n");
    spamm_print_dense(N[0], N[1], row_major, B_dense);
#endif

    /* Convert to SpAMM matrix. */
    A = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, A_dense);
    B = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, B_dense);

    if(spamm_check(A, TEST_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }
    if(spamm_check(B, TEST_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        if(A_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
        {
          result |= SPAMM_ERROR;
          printf("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
              spamm_get(i, A), A_dense[i[0]*N[1]+i[1]]);
          break;
        }

        if(B_dense[i[0]*N[1]+i[1]] != spamm_get(i, B))
        {
          result |= SPAMM_ERROR;
          printf("failed test at B[%u][%u] (found %f, should be %f)\n", i[0], i[1],
              spamm_get(i, B), B_dense[i[0]*N[1]+i[1]]);
          break;
        }
      }
    }

    /* Add by hand for verification. */
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] =
          alpha*A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]
          + beta*B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])];
      }
    }

#ifdef PRINT_MATRIX
    printf("A_dense =\n");
    spamm_print_dense(N[0], N[1], row_major, A_dense);
#endif

    flop = 0;
    spamm_add(alpha, A, beta, B, &flop);
    printf("%e flop\n", flop);

#ifdef PRINT_MATRIX
    printf("A =\n");
    spamm_print(A);
#endif

    /* Check tree consistency. */
    if(spamm_check(A, TEST_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    /* Compare result. */
    max_diff = 0.0;
    max_diff_i[0] = 0;
    max_diff_i[1] = 0;
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        if(fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A)) > max_diff)
        {
          max_diff = fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A));
          max_diff_i[0] = i[0];
          max_diff_i[1] = i[1];
        }
      }
    }

    /* Test result. */
    if(max_diff > TEST_TOLERANCE)
    {
      result |= SPAMM_ERROR;
      printf("max diff of A[%u][%u] (sparse_A = %u, sparse_B = %u, test tolerance was %1.2e) = %e\n",
          max_diff_i[0], max_diff_i[1], sparse_A[sparse_test], sparse_B[sparse_test], TEST_TOLERANCE, max_diff);
    }

    else
    {
      printf("passed\n");
    }

    free(A_dense);
    free(B_dense);

    spamm_delete(&A);
    spamm_delete(&B);
  }

  return result;
}

#include "spamm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_TOLERANCE 5e-7

//#define PRINT_MATRIX

void
fill_matrix (const unsigned int number_dimensions,
    const short is_sparse,
    const unsigned int *const N,
    float *A_dense)
{
  int dim;
  unsigned int *i;
  unsigned int N_contiguous;

  for (dim = 0, N_contiguous = 1; dim < number_dimensions; dim++)
  {
    N_contiguous *= N[dim];
  }

  i = calloc(number_dimensions, sizeof(unsigned int));

  if (is_sparse)
  {
    switch (number_dimensions)
    {
      case 1:
        for (i[0] = 0; i[0] < N[0]; i[0]++)
        {
          A_dense[i[0]] = rand()/(float) RAND_MAX;
        }
        break;

      case 2:
        for (i[0] = 0; i[0] < N[0]; i[0]++) {
          for ((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
          {
            A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(float) RAND_MAX;
          }
        }
        break;

      case 3:
        for (i[0] = 0; i[0] < N[0]; i[0]++) {
          for ((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++) {
            for ((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < N[2]; i[2]++)
            {
              A_dense[i[0]+N[0]*(i[1]+N[1]*i[2])] = rand()/(float) RAND_MAX;
            }
          }
        }
        break;

      default:
        SPAMM_FATAL("FIXME\n");
        break;
    }
  }

  else
  {
    for (i[0] = 0; i[0] < N_contiguous; i[0]++)
    {
      A_dense[i[0]] = rand()/(float) RAND_MAX;
    }
  }

  free(i);
}

int
main (int argc, char **argv)
{
  int result = 0;

  unsigned int number_dimensions;

  unsigned int *N;
  unsigned int *i;

  unsigned int N_contiguous;

  unsigned int dim;

  const unsigned int chunk_tier = 4;

  short use_linear_tree;
  short is_sparse_A;
  short is_sparse_B;

  double alpha = 1.2;
  double beta = 0.8;

  float *A_dense;
  float *B_dense;

  float max_diff;
  unsigned int *max_diff_i;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;

  for (number_dimensions = 1; number_dimensions <= 3; number_dimensions++)
  {
    i = calloc(number_dimensions, sizeof(unsigned int));
    N = calloc(number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < number_dimensions; dim++)
    {
      N[dim] = 150+(int) ((0.5-(float) rand()/(float) RAND_MAX)*30);
    }

    N_contiguous = 1;
    for (dim = 0; dim < number_dimensions; dim++)
    {
      N_contiguous *= N[dim];
    }

    A_dense = (float*) calloc(N_contiguous, sizeof(float));
    B_dense = (float*) calloc(N_contiguous, sizeof(float));

    for (use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      for (is_sparse_A = 0; is_sparse_A < 2; is_sparse_A++) {
        for (is_sparse_B = 0; is_sparse_B < 2; is_sparse_B++)
        {
          printf("dim: %u, linTree: %u, sparse A: %u, sparse B: %u\n",
              number_dimensions, use_linear_tree, is_sparse_A, is_sparse_B);

          fill_matrix(number_dimensions, is_sparse_A, N, A_dense);
          fill_matrix(number_dimensions, is_sparse_B, N, B_dense);

          /* Convert to SpAMM matrix. */
          A = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, A_dense);
          B = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, B_dense);

          printf("A: ");
          spamm_print_info(A);

          printf("B: ");
          spamm_print_info(B);

          /* Add by hand for verification. */
          switch (number_dimensions)
          {
            case 1:
              for (i[0] = 0; i[0] < N[0]; i[0]++)
              {
                A_dense[i[0]] = alpha*A_dense[i[0]] + beta*B_dense[i[0]];
              }
              break;

            case 2:
              for (i[0] = 0; i[0] < N[0]; i[0]++) {
                for (i[1] = 0; i[1] < N[1]; i[1]++)
                {
                  A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] =
                    alpha*A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]
                    + beta*B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])];
                }
              }
              break;

            case 3:
              for (i[0] = 0; i[0] < N[0]; i[0]++) {
                for (i[1] = 0; i[1] < N[1]; i[1]++) {
                  for (i[2] = 0; i[2] < N[2]; i[2]++)
                  {
                    A_dense[spamm_index_row_major_3(number_dimensions, N, i)] =
                      alpha*A_dense[spamm_index_row_major_3(number_dimensions, N, i)]
                      + beta*B_dense[spamm_index_row_major_3(number_dimensions, N, i)];
                  }
                }
              }
              break;

            default:
              SPAMM_FATAL("FIXME\n");
              break;
          }

          spamm_add(alpha, A, beta, B);

          /* Check tree consistency. */
          result |= spamm_check(A, TEST_TOLERANCE);

          /* Compare result. */
          max_diff = 0.0;
          max_diff_i = calloc(number_dimensions, sizeof(unsigned int));

          switch (number_dimensions)
          {
            case 1:
              for (i[0] = 0; i[0] < N[0]; i[0]++)
              {
                if (fabs(A_dense[i[0]]-spamm_get(i, A)) > max_diff)
                {
                  max_diff = fabs(A_dense[i[0]]-spamm_get(i, A));
                  max_diff_i[0] = i[0];
                }
              }

              if (max_diff > TEST_TOLERANCE)
              {
                result |= SPAMM_ERROR;
                printf("[add_spamm] max diff of A[%u] (test tolerance was %1.2e) = %e\n",
                    max_diff_i[0], TEST_TOLERANCE, max_diff);
              }
              break;

            case 2:
              for (i[0] = 0; i[0] < N[0]; i[0]++) {
                for (i[1] = 0; i[1] < N[1]; i[1]++)
                {
                  if (fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A)) > max_diff)
                  {
                    max_diff = fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A));
                    max_diff_i[0] = i[0];
                    max_diff_i[1] = i[1];
                  }
                }
              }

              if (max_diff > TEST_TOLERANCE)
              {
                result |= SPAMM_ERROR;
                printf("[add_spamm] max diff of A[%u][%u] (test tolerance was %1.2e) = %e\n",
                    max_diff_i[0], max_diff_i[1], TEST_TOLERANCE, max_diff);
              }
              break;

            case 3:
              for (i[0] = 0; i[0] < N[0]; i[0]++) {
                for (i[1] = 0; i[1] < N[1]; i[1]++) {
                  for (i[2] = 0; i[2] < N[2]; i[2]++)
                  {
                    if (fabs(A_dense[spamm_index_row_major_3(number_dimensions, N, i)]-spamm_get(i, A)) > max_diff)
                    {
                      max_diff = fabs(A_dense[spamm_index_row_major_3(number_dimensions, N, i)]-spamm_get(i, A));
                      max_diff_i[0] = i[0];
                      max_diff_i[1] = i[1];
                      max_diff_i[2] = i[2];
                    }
                  }
                }
              }

              if (max_diff > TEST_TOLERANCE)
              {
                result |= SPAMM_ERROR;
                printf("[add_spamm] max diff of A[%u][%u][%u] (test tolerance was %1.2e) = %e\n",
                    max_diff_i[0], max_diff_i[1], max_diff_i[2], TEST_TOLERANCE, max_diff);
              }
              break;

            default:
              SPAMM_FATAL("FIXME\n");
              break;
          }
          free(max_diff_i);
          spamm_delete(&A);
          spamm_delete(&B);
        }
      }
    }
    free(i);
    free(N);
    free(A_dense);
    free(B_dense);
  }

  return result;
}

#include "spamm.h"

#include "test.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REL_TOLERANCE 1e-8

#define TEST_TOLERANCE 5e-7

//#define PRINT_MATRIX

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

  double flop;

  for(number_dimensions = 1; number_dimensions <= 3; number_dimensions++) {
    for(use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++) {
      for(is_sparse_A = 0; is_sparse_A < 2; is_sparse_A++) {
        for(is_sparse_B = 0; is_sparse_B < 2; is_sparse_B++)
        {
          printf("dim: %u, linTree: %u, sparse A: %u, sparse B: %u\n",
              number_dimensions, use_linear_tree, is_sparse_A, is_sparse_B);

          i = calloc(number_dimensions, sizeof(unsigned int));
          N = generate_shape(number_dimensions, 0);
          N_contiguous = 1;
          for(dim = 0; dim < number_dimensions; dim++)
          {
            N_contiguous *= N[dim];
          }

          A_dense = generate_matrix(number_dimensions, is_sparse_A, N);
          B_dense = generate_matrix(number_dimensions, is_sparse_B, N);

          /* Convert to SpAMM matrix. */
          A = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, A_dense);
          B = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, B_dense);

          printf("A: ");
          spamm_print_info(A);

          printf("B: ");
          spamm_print_info(B);

          /* Add by hand for verification. */
          switch(number_dimensions)
          {
            case 1:
              for(i[0] = 0; i[0] < N[0]; i[0]++)
              {
                A_dense[i[0]] = alpha*A_dense[i[0]] + beta*B_dense[i[0]];
              }
              break;

            case 2:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for(i[1] = 0; i[1] < N[1]; i[1]++)
                {
                  A_dense[spamm_index_row_major_3(number_dimensions, N, i)] =
                    alpha*A_dense[spamm_index_row_major_3(number_dimensions, N, i)]
                    + beta*B_dense[spamm_index_row_major_3(number_dimensions, N, i)];
                }
              }
              break;

            case 3:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for(i[1] = 0; i[1] < N[1]; i[1]++) {
                  for(i[2] = 0; i[2] < N[2]; i[2]++)
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

          flop = 0;
          spamm_add(alpha, A, beta, B, &flop);
          printf("[add_spamm] %e flop\n", flop);

          /* Check tree consistency. */
          if(spamm_check(A, REL_TOLERANCE) != SPAMM_OK)
          {
            SPAMM_FATAL("failed\n");
          }

          if(compare_spamm_to_dense(A, A_dense, TEST_TOLERANCE) != SPAMM_OK)
          {
            SPAMM_FATAL("comparison failed\n");
          }

          free(i);
          free(N);
          free(A_dense);
          free(B_dense);
          free(max_diff_i);
          spamm_delete(&A);
          spamm_delete(&B);
        }
      }
    }
  }

  return result;
}

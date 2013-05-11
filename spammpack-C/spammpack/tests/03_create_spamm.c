#include "test.h"

#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

#define REL_TOLERANCE 1e-8
#define ABS_TOLERANCE 1e-8
#define DECAY 2.5

int
main (int argc, char **argv)
{
  int result = SPAMM_OK;

  unsigned int number_dimensions;

  unsigned int *N;

  const unsigned int N_mean = 150;
  const unsigned int N_width = 30;

  const unsigned int chunk_tier = 3;

  short use_linear_tree;

  enum matrix_t matrix_type;

  struct spamm_matrix_t *A;

  float *A_dense;

  for(number_dimensions = 1; number_dimensions <= 3; number_dimensions++) {
    for(use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      if(number_dimensions != 2 && use_linear_tree)
      {
        continue;
      }

      for(matrix_type = 0; matrix_type < NUMBER_MATRIX_TYPES; matrix_type++)
      {
        printf("dim: %u, linTree: %u, matrix_type: %s\n", number_dimensions, use_linear_tree, get_matrix_type_name(matrix_type));

        if(matrix_type == exponential_decay)
        {
          N = generate_shape(number_dimensions, 1, N_mean, N_width);
        }

        else
        {
          N = generate_shape(number_dimensions, 0, N_mean, N_width);
        }

        A_dense = generate_matrix_float(number_dimensions, matrix_type, N, DECAY);
        A = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, A_dense);

        //printf("A info: ");
        //spamm_print_info(A);

        if(spamm_check(A, REL_TOLERANCE) != SPAMM_OK)
        {
          SPAMM_FATAL("failed\n");
        }

        if(compare_spamm_to_dense_float(A, A_dense, ABS_TOLERANCE) != SPAMM_OK)
        {
          SPAMM_FATAL("comparison failed\n");
        }

#ifdef PRINT_DEBUG
        SPAMM_INFO("test passed\n");
#endif

        free(A_dense);
        free(N);
        spamm_delete(&A);
      }
    }
  }

  return result;
}

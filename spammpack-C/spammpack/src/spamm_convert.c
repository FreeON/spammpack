/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Convert a dense matrix to SpAMM.
 *
 * @param number_dimensions The number of dimensions.
 * @param N The number of rows/columns.
 * @param chunk_tier The tier at which to store contiguous submatrix
 * blocks in the hierarhical tree.
 * @param use_linear_tree If set to zero, then the tree will be stored in the
 * hierachical format, otherwise storage will switch to linear format at
 * chunk_tier.
 * @param dense_type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 * @param spamm_layout The layout of the SpAMM data nodes.
 *
 * @return The SpAMM matrix.
 */
struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const enum spamm_layout_t dense_type,
    const float *const A_dense)
{
  struct spamm_matrix_t *A;
  unsigned int *i;

  assert(A_dense != NULL);

  A = spamm_new(number_dimensions, N, chunk_tier, use_linear_tree);

  i = calloc(number_dimensions, sizeof(unsigned int));

  switch(number_dimensions)
  {
    case 1:
      for(i[0] = 0; i[0] < N[0]; i[0]++)
      {
        switch(dense_type)
        {
          case row_major:
          case column_major:
            spamm_set(i, A_dense[i[0]], A);
            break;

          default:
            SPAMM_FATAL("unknown type\n");
            break;
        }
      }
      break;

    case 2:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++)
        {
          switch(dense_type)
          {
            case row_major:
              spamm_set(i, A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])], A);
              break;

            case column_major:
              spamm_set(i, A_dense[spamm_index_column_major(i[0], i[1], N[0], N[1])], A);
              break;

            default:
              SPAMM_FATAL("unknown type\n");
              break;
          }
        }
      }
      break;

    case 3:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++) {
          for(i[2] = 0; i[2] < N[2]; i[2]++)
          {
            switch(dense_type)
            {
              case row_major:
                spamm_set(i, A_dense[spamm_index_row_major_3(number_dimensions, N, i)], A);
                break;

              case column_major:
                spamm_set(i, A_dense[spamm_index_column_major_3(number_dimensions, N, i)], A);
                break;

              default:
                SPAMM_FATAL("unknown type\n");
                break;
            }
          }
        }
      }
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }
  free(i);

  return A;
}

/** Convert a SpAMM matrix to a dense matrix.
 *
 * @param A The SpAMM matrix.
 *
 * @return The newly allocated dense matrix. This matrix has to be freed by
 * the caller when it is not needed anymore.
 */
float *
spamm_convert_spamm_to_dense (const struct spamm_matrix_t *const A)
{
  float *ADense = NULL;
  unsigned int *i;

  switch(A->number_dimensions)
  {
    case 2:
      ADense = calloc(A->N[0]*A->N[1], sizeof(float));
      i = calloc(2, sizeof(unsigned int));
      for (i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for (i[1] = 0; i[1] < A->N[1]; i[1]++)
        {
          ADense[spamm_index_column_major_3(A->number_dimensions, A->N, i)] = spamm_get(i, A);
        }
      }
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }

  return ADense;
}

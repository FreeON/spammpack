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
 * @param contiguous_tier The tier at which to store contiguous submatrix
 * blocks in the hierarhical tree.
 * @param N_block The size of matrix to which the SpAMM condition is applied.
 * @param use_linear_tree If set to zero, then the tree will be stored in the
 * hierachical format, otherwise storage will switch to linear format at
 * contiguous_tier.
 * @param dense_type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 * @param spamm_layout The layout of the SpAMM data nodes.
 *
 * @return The SpAMM matrix.
 */
struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int contiguous_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const enum spamm_layout_t dense_type,
    const float *const A_dense,
    const enum spamm_layout_t spamm_layout)
{
  struct spamm_matrix_t *A = NULL;
  unsigned int *i;

  assert(A_dense != NULL);

  if (number_dimensions != 2)
  {
    SPAMM_FATAL("can not handle this case\n");
  }

  A = spamm_new(number_dimensions, N, contiguous_tier, N_block, use_linear_tree, spamm_layout);

  i = calloc(number_dimensions, sizeof(unsigned int));
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
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
  free(i);

  return A;
}

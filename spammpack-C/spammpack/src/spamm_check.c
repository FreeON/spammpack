#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float
spamm_chunk_check (spamm_chunk_t *chunk,
    const float tolerance)
{
  float norm2_reference;

  float *norm;
  float *norm2;
  float *matrix;

  unsigned i, j;
  unsigned int tier;
  unsigned int N_contiguous;
  unsigned int N_block;
  unsigned int use_linear_tree;

  if (chunk == NULL) { return 0.0; }

  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);

  matrix = spamm_chunk_get_matrix(chunk);
  norm = spamm_chunk_get_norm(chunk);
  norm2 = spamm_chunk_get_norm2(chunk);

  norm2_reference = 0;
  for (i = 0; i < N_contiguous; i++) {
    for (j = 0; j < N_contiguous; j++)
    {
      norm2_reference += matrix[i*N_contiguous+j]*matrix[i*N_contiguous+j];
    }
  }

  if (fabs(norm2_reference - norm2[0]) > tolerance)
  {
    SPAMM_WARN("norm2 mismatch: found %e, should be %e (abs. diff = %e)\n",
        norm2[0], norm2_reference, fabs(norm2[0]-norm2_reference));
  }

  if (fabs(sqrt(norm2_reference) - norm[0]) > tolerance)
  {
    SPAMM_WARN("norm mismatch\n");
  }

  if (use_linear_tree)
  {
  }

  return norm2[0];
}

float
spamm_recursive_check (const struct spamm_recursive_node_t *const node,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const float tolerance)
{
  int child_index;

  float norm2;

  if (node == NULL) { return SPAMM_OK; }

  if (tier == chunk_tier)
  {
    norm2 = spamm_chunk_check(node->tree.chunk, tolerance);
  }

  else
  {
    norm2 = 0;
    for (child_index = 0; child_index < ipow(2, number_dimensions); child_index++)
    {
      norm2 += spamm_recursive_check(node->tree.child[child_index], number_dimensions, tier+1, chunk_tier, tolerance);
    }
  }

  if (fabs(norm2 - node->norm2) > tolerance)
  {
    SPAMM_WARN("norm2 mismatch\n");
  }

  if (fabs(sqrt(norm2) - node->norm) > tolerance)
  {
    SPAMM_WARN("norm mismatch\n");
  }

  return norm2;
}

/** Check the internal consistency of a matrix.
 *
 * @param A The matrix to check
 * @param tolerance The absolute tolerance when comparing values.
 *
 * @return The following error codes are returned:
 *   - SPAMM_OK - The matrix is consistent.
 *   - SPAMM_ERROR - Something is not consistent.
 */
void
spamm_check (const struct spamm_matrix_t *A, const float tolerance)
{
  float norm2;

  if (A->chunk_tier == 0)
  {
    norm2 = spamm_chunk_check(A->tree.chunk, tolerance);
  }

  else
  {
    norm2 = spamm_recursive_check(A->tree.recursive_tree, A->number_dimensions, 0, A->chunk_tier, tolerance);
  }
}

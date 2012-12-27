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

  unsigned i, j, k;
  unsigned int tier;
  unsigned int number_dimensions;
  unsigned int N_contiguous;
  unsigned int use_linear_tree;
  unsigned int number_tiers;
  unsigned int linear_index;
  unsigned int offset;

  if(chunk == NULL) { return 0.0; }

  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  matrix = spamm_chunk_get_matrix(chunk);

  if(use_linear_tree)
  {
    for(tier = 0; tier < number_tiers-SPAMM_KERNEL_DEPTH; tier++)
    {
      norm = spamm_chunk_get_tier_norm(tier, chunk);
      norm2 = spamm_chunk_get_tier_norm2(tier, chunk);

      for(linear_index = 0; linear_index < ipow(4, tier); linear_index++)
      {
        for(i = linear_index*ipow(N_contiguous >> tier, 2), norm2_reference = 0;
            i < (linear_index+1)*ipow(N_contiguous >> tier, 2); i++)
        {
          norm2_reference += matrix[i]*matrix[i];
        }

        if(fabs(norm2_reference - norm2[linear_index]) > tolerance)
        {
          SPAMM_WARN("norm2 mismatch: tier = %u, linear index = %u, found %e, should be %e (abs. diff = %e, rel. diff = %e))\n",
              tier, linear_index,
              norm2[linear_index], norm2_reference,
              fabs(norm2[linear_index]-norm2_reference),
              fabs(norm2[linear_index]-norm2_reference)/norm2_reference);
        }

        if(fabs(sqrt(norm2_reference) - norm[linear_index]) > tolerance)
        {
          SPAMM_WARN("norm mismatch: tier = %u, linear index = %u, found %e, should be %e (abs. diff = %e, rel. diff = %e))\n",
              tier, linear_index,
              norm[linear_index], sqrt(norm2_reference),
              fabs(norm[linear_index]-sqrt(norm2_reference)),
              fabs(norm[linear_index]-sqrt(norm2_reference))/sqrt(norm2_reference));
        }
      }
    }

    /* Check last tier (SPAMM_N_BLOCK). */
    norm = spamm_chunk_get_tier_norm(number_tiers-1, chunk);
    norm2 = spamm_chunk_get_tier_norm2(number_tiers-1, chunk);

    for(linear_index = 0, norm2_reference = 0;
        linear_index < ipow(4, number_tiers-number_tiers-1);
        linear_index++)
    {
      for(i = 0; i < SPAMM_N_KERNEL_BLOCKED; i++) {
        for(j = 0; j < SPAMM_N_KERNEL_BLOCKED; j++)
        {
          offset = linear_index*ipow(SPAMM_N_KERNEL, 2)
            + (i*SPAMM_N_KERNEL_BLOCKED+j)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;
          for(k = 0; k < SPAMM_N_BLOCK*SPAMM_N_BLOCK; k++)
          {
            norm2_reference += matrix[offset+k]*matrix[offset+k];
          }

          if(fabs(norm2_reference - norm2[linear_index]) > tolerance)
          {
            SPAMM_WARN("norm2 mismatch: tier = %u, linear_index = %u, found %e, should be %e (abs. diff = %e, rel. diff = %e)\n",
                number_tiers-1, linear_index,
                norm2[linear_index], norm2_reference,
                fabs(norm2[linear_index]-norm2_reference),
                fabs(norm2[linear_index]-norm2_reference)/norm2_reference);
          }

          if(fabs(sqrt(norm2_reference) - norm[linear_index]) > tolerance)
          {
            SPAMM_WARN("norm mismatch\n");
          }
        }
      }
    }
  }

  norm = spamm_chunk_get_norm(chunk);
  norm2 = spamm_chunk_get_norm2(chunk);

  norm2_reference = 0;
  for(i = 0; i < ipow(N_contiguous, number_dimensions); i++)
  {
    norm2_reference += matrix[i]*matrix[i];
  }

  if(fabs(norm2_reference - norm2[0]) > tolerance)
  {
    SPAMM_WARN("norm2 mismatch: found %e, should be %e (abs. diff = %e, rel. diff = %e)\n",
        norm2[0],
        norm2_reference,
        fabs(norm2[0]-norm2_reference),
        fabs(norm2[0]-norm2_reference)/norm2_reference);
  }

  if(fabs(sqrt(norm2_reference) - norm[0]) > tolerance)
  {
    SPAMM_WARN("norm mismatch: found %e, should be %e (abs. diff = %e, rel. diff = %e)\n",
        norm[0],
        sqrt(norm2_reference),
        fabs(norm[0]-sqrt(norm2_reference)),
        fabs(norm[0]-sqrt(norm2_reference))/sqrt(norm2_reference));
  }

  return norm2_reference;
}

float
spamm_recursive_check (const struct spamm_recursive_node_t *const node,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const float tolerance)
{
  int child_index;

  float norm2_reference;

  if(node == NULL) { return SPAMM_OK; }

  if(tier == chunk_tier)
  {
    norm2_reference = spamm_chunk_check(node->tree.chunk, tolerance);
  }

  else
  {
    norm2_reference = 0;
    for(child_index = 0; child_index < ipow(2, number_dimensions); child_index++)
    {
      norm2_reference += spamm_recursive_check(node->tree.child[child_index],
          number_dimensions, tier+1, chunk_tier, tolerance);
    }
  }

  if(fabs(norm2_reference - node->norm2) > tolerance)
  {
    SPAMM_WARN("norm2 mismatch: found %e, should be %e (abs. diff = %e, rel. diff = %e)\n",
        node->norm2,
        norm2_reference,
        fabs(node->norm2-norm2_reference),
        fabs(node->norm2-norm2_reference)/norm2_reference);
  }

  if(fabs(sqrt(norm2_reference) - node->norm) > tolerance)
  {
    SPAMM_WARN("norm mismatch: found %e, should be %e (abs. diff = %e, rel. diff = %e)\n",
        node->norm,
        sqrt(norm2_reference),
        fabs(node->norm-sqrt(norm2_reference)),
        fabs(node->norm-sqrt(norm2_reference))/sqrt(norm2_reference));
  }

  return norm2_reference;
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

  norm2 = spamm_recursive_check(A->recursive_tree, A->number_dimensions, 0, A->chunk_tier, tolerance);
}

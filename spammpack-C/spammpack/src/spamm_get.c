#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Get an element from a SpAMM chunk.
 *
 * @param i The row/column index array.
 * @param chunk The chunk.
 *
 * @return The matrix element.
 */
float
spamm_chunk_get (const unsigned int *i,
    spamm_chunk_t *chunk)
{
  float Aij = 0;

  int dim;

  short use_linear_tree;

  unsigned int tier;
  unsigned int linear_index;

  unsigned int number_dimensions;
  unsigned int number_tiers;

  unsigned int *N_lower;
  unsigned int *N_upper;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  unsigned int *new_i;

  float *norm;
  float *norm2;

  float *A;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);

  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);

  new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
  new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    new_N_lower[dim] = N_lower[dim];
    new_N_upper[dim] = N_upper[dim];
  }

  /* Correct tier count. */
  if (use_linear_tree)
  {
    number_tiers -= SPAMM_KERNEL_DEPTH;
  }

  /* Z-curve ordering down to SPAMM_N_KERNEL. */
  for (tier = 0, linear_index = 0; tier < number_tiers; tier++)
  {
    norm = spamm_chunk_get_tier_norm(tier, chunk);
    norm2 = spamm_chunk_get_tier_norm2(tier, chunk);

    /* Update norm. */
    norm2[linear_index] += Aij*Aij;
    norm[linear_index] = sqrt(norm2[linear_index]);

    if (tier+1 < number_tiers)
    {
      /* Recurse. */
      linear_index <<= number_dimensions;

      for (dim = 0; dim < number_dimensions; dim++)
      {
        if (i[dim] < new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2)
        {
          new_N_lower[dim] = new_N_lower[dim];
          new_N_upper[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
        }

        else
        {
          new_N_upper[dim] = new_N_upper[dim];
          new_N_lower[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
          linear_index |= (1 << dim);
        }
      }
    }

    else
    {
      break;
    }
  }

  A = spamm_chunk_get_matrix(chunk);
  if (use_linear_tree)
  {
    Aij = A[linear_index*ipow(SPAMM_N_KERNEL, number_dimensions)*sizeof(float)
      +spamm_index_kernel_block(i[0]-new_N_lower[0], i[1]-new_N_lower[1], row_major)];
  }

  else
  {
    new_i = calloc(number_dimensions, sizeof(unsigned int));
    for (dim = 0; dim < number_dimensions; dim++)
    {
      new_i[dim] = i[dim]-N_lower[dim];
    }
    Aij = A[spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], new_i)];
    free(new_i);
  }

  free(new_N_lower);
  free(new_N_upper);

  return Aij;
}

/** Get an element from a recursive matrix.
 *
 * If the matrix is NULL, this function returns 0.
 *
 * @param i An array of row indices.
 * @param A The matrix.
 *
 * @return The matrix element Aij.
 */
float
spamm_recursive_get (const unsigned int number_dimensions,
    const unsigned int *const i,
    const unsigned int *N_lower,
    const unsigned int *N_upper,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const struct spamm_recursive_node_t *node)
{
  int dim;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  short child_index;

  float Aij;

  if (node == NULL) { return 0; }

  if (tier == chunk_tier)
  {
    Aij = spamm_chunk_get(i, node->tree.chunk);
  }

  else
  {
    new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

    child_index = 0;

    for (dim = 0; dim < number_dimensions; dim++)
    {
      if (i[dim] < N_lower[dim]+(N_upper[dim]-N_lower[dim])/2)
      {
        new_N_lower[dim] = N_lower[dim];
        new_N_upper[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
      }

      else
      {
        new_N_lower[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
        new_N_upper[dim] = N_upper[dim];
        child_index |= (1 << dim);
      }
    }

    Aij = spamm_recursive_get(number_dimensions, i, new_N_lower, new_N_upper,
        tier+1, chunk_tier, use_linear_tree, node->tree.child[child_index]);

    free(new_N_lower);
    free(new_N_upper);
  }

  return Aij;
}

/** Get an element from a matrix.
 *
 * @param i The row/column index.
 * @param A The matrix.
 *
 * @return The matrix element.
 */
float
spamm_get (const unsigned int *const i, const struct spamm_matrix_t *A)
{
  unsigned int *N_lower;
  unsigned int *N_upper;

  int dim;

  float Aij;

  if (A->chunk_tier == 0)
  {
    Aij = spamm_chunk_get(i, A->tree.chunk);
  }

  else
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    Aij = spamm_recursive_get(A->number_dimensions, i, N_lower, N_upper, 0,
        A->chunk_tier, A->use_linear_tree, A->tree.recursive_tree);

    free(N_lower);
    free(N_upper);
  }

  return Aij;
}

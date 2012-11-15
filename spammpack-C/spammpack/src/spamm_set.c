#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW_NORM
#define SPAMM_SET_NO_ZERO

/** Set an element in a SpAMM chunk.
 *
 * @param i The row/column index array.
 * @param Aij The value of the matrix element.
 * @param chunk The SpAMM chunk.
 */
void
spamm_chunk_set (const unsigned int *const i,
    const float Aij,
    spamm_chunk_t *chunk)
{
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
    A[linear_index*ipow(SPAMM_N_KERNEL, number_dimensions)*sizeof(float)
      +spamm_index_kernel_block(i[0]-new_N_lower[0], i[1]-new_N_lower[1], row_major)] = Aij;
  }

  else
  {
    new_i = calloc(number_dimensions, sizeof(unsigned int));
    for (dim = 0; dim < number_dimensions; dim++)
    {
      new_i[dim] = i[dim]-N_lower[dim];
    }
    A[spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], new_i)] = Aij;
    free(new_i);
  }

  free(new_N_lower);
  free(new_N_upper);
}

/** Recursively set a matrix element.
 *
 * @param number_dimensions The number of dimensions.
 * @param i An array of row/column indices.
 * @param Aij The value of the matrix element A(i,j).
 * @param N An array of matrix dimensions (unpadded).
 * @param N_lower An array of left-most column indices.
 * @param N_upper An array of right-most column indices.
 * @param tier The tier this node is on.
 * @param chunk_tier The size of the contiguous submatrix block.
 * @param linear_tier The size of the submatrix that is stored in hashed format.
 * @param layout The layout of the matrix elements.
 * @param node The node.
 */
void
spamm_recursive_set (const unsigned int number_dimensions,
    const unsigned int *const i,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const unsigned int depth,
    const float Aij,
    struct spamm_recursive_node_t **node)
{
  int dim;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  short child_index;

  if (*node == NULL)
  {
    *node = spamm_recursive_new_node();
  }

  /* Update norm. */
  (*node)->norm2 += Aij*Aij;
  (*node)->norm   = sqrt((*node)->norm2);

  if (tier == chunk_tier)
  {
    if ((*node)->tree.chunk == NULL)
    {
      (*node)->tree.chunk = spamm_new_chunk(number_dimensions,
          use_linear_tree, N, N_lower, N_upper);
    }

    spamm_chunk_set(i, Aij, (*node)->tree.chunk);
  }

  else
  {
    if ((*node)->tree.child == NULL)
    {
      (*node)->tree.child = calloc(ipow(2, number_dimensions), sizeof(struct spamm_recursive_node_t*));
    }

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

    spamm_recursive_set(number_dimensions, i, N, new_N_lower, new_N_upper,
        tier+1, chunk_tier, use_linear_tree, depth, Aij,
        &(*node)->tree.child[child_index]);

    free(new_N_lower);
    free(new_N_upper);
  }
}

/** Set an element in a matrix.
 *
 * @param i The row/column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param A The matrix.
 */
void
spamm_set (const unsigned int *const i, const float Aij, struct spamm_matrix_t *A)
{
  unsigned int *N_lower;
  unsigned int *N_upper;

  int dim;

  if (A->chunk_tier == 0)
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    if (A->tree.chunk == NULL)
    {
      A->tree.chunk = spamm_new_chunk(A->number_dimensions,
          A->use_linear_tree, A->N, N_lower, N_upper);
    }

    spamm_chunk_set(i, Aij, A->tree.chunk);

    free(N_lower);
    free(N_upper);
  }

  else
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    spamm_recursive_set(A->number_dimensions, i, A->N, N_lower, N_upper, 0,
        A->chunk_tier, A->use_linear_tree, A->depth, Aij,
        &(A->tree.recursive_tree));

    free(N_lower);
    free(N_upper);
  }
}

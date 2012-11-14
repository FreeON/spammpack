#include "spamm.h"
#include "spamm_types_private.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Get the depth of a matrix tree.
 *
 * @return The depth.
 */
unsigned int
spamm_get_tree_depth (const unsigned int number_dimensions,
    const unsigned int *const N,
    const short use_linear_tree)
{
  int dim;
  unsigned int depth;
  unsigned int N_temp;
  double x, x_N;

  /* Pad to powers of M_child x N_child. */
  x = 0;
  for (dim = 0; dim < number_dimensions; dim++)
  {
    /* Make sure we pad at least to the extend that we can store that linear
     * kernel matrix. */
    if (use_linear_tree)
    {
      if (N[dim] < SPAMM_N_KERNEL)
      {
        N_temp = SPAMM_N_KERNEL;
      }

      else
      {
        N_temp = N[dim];
      }
    }

    else
    {
      N_temp = N[dim];
    }

    x_N = log(N_temp)/log(2);
    if (x_N > x)
    {
      x = x_N;
    }
  }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (depth >= 1)
  {
    for (dim = 0; dim < number_dimensions; dim++)
    {
      if ((int) (ipow(2, depth-1)) < N[dim])
      {
        depth++;
        break;
      }
    }
    depth--;
  }

  return depth;
}

/** Allocate a new node of a recursive matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param number_dimensions The number of dimensions.
 * @param chunk_tier The tier at which to store contiguous submatrix
 * blocks in the hierarhical tree.
 * @param N_block The size of matrix to which the SpAMM condition is applied.
 * @param use_linear_tree If set to zero, then the tree will be stored in the
 * hierachical format, otherwise storage will switch to linear format at
 * chunk_tier.
 * @param N An array of the matrix dimensions (unpadded).
 * @param N_lower An array of the lowest row index of this submatrix node.
 * @param N_upper An array of the lowest row index of this submatrix node.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_recursive_node_t *
spamm_recursive_new_node ()
{
  struct spamm_recursive_node_t *node = NULL;

  /* Allocate memory. */
  node = calloc(1, sizeof(struct spamm_recursive_node_t));

  return node;
}

/** Initialize a new matrix object.
 *
 * Two settings determine when and if the tree is stored in linear or
 * hierarchical format. If N_linear == N_contiguous, then the tree format will
 * switch to linear tree format at N_contiguous. If N_linear < N_contiguous,
 * then the whole tree will be stored hierarchically and the spamm condition
 * is applied to the SpAMM chunks.
 *
 * @param number_dimensions The number of dimensions of this matrix.
 * @param N The number of rows/columns of the matrix. This array has to have
 * a size of number_dimensions.
 * @param chunk_tier The tier at which to store contiguous submatrix
 * blocks in the hierarhical tree.
 * @param use_linear_tree If set to zero, then the tree will be stored in the
 * hierachical format, otherwise storage will switch to linear format at
 * chunk_tier.
 *
 * @return The newly allocated matrix. This matrix has to be freed by calling
 * spamm_delete().
 */
struct spamm_matrix_t *
spamm_new (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree)
{
  int dim;
  struct spamm_matrix_t *A = NULL;

  for (dim = 0; dim < number_dimensions; dim++)
  {
    if (N[dim] == 0)
    {
      SPAMM_FATAL("N[%u] == 0\n", dim);
    }
  }

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_matrix_t));

  /* Store the number of dimensions. */
  A->number_dimensions = number_dimensions;

  /* Store matrix dimensions. */
  A->N = calloc(number_dimensions, sizeof(unsigned int));
  for (dim = 0; dim < number_dimensions; dim++)
  {
    A->N[dim] = N[dim];
  }

  if (number_dimensions == 2 && use_linear_tree)
  {
    A->use_linear_tree = use_linear_tree;
  }

  /* Get tree depth. */
  A->depth = spamm_get_tree_depth(number_dimensions, A->N, use_linear_tree);

  /* Set padded matrix size. */
  A->N_padded = (int) (ipow(2, A->depth));

  /* Adjust the depth. */
  if (number_dimensions == 2 && use_linear_tree)
  {
    A->depth -= 3; /* 16x16 submatrix blocks for linear kernel. */
  }

  if (chunk_tier > A->depth)
  {
    SPAMM_WARN("chunk tier (%u) is greater than depth (%u)\n", chunk_tier, A->depth);
    A->chunk_tier = A->depth;
  }

  else
  {
    A->chunk_tier = chunk_tier;
    A->depth = chunk_tier;
  }

  /* Done. */
  return A;
}

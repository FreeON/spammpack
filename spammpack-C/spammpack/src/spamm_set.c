#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW_NORM
#define SPAMM_SET_NO_ZERO

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
    const unsigned int kernel_tier,
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
        tier+1, chunk_tier, kernel_tier, use_linear_tree, depth, Aij,
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
        A->chunk_tier, A->kernel_tier, A->use_linear_tree, A->depth, Aij,
        &(A->tree.recursive_tree));

    free(N_lower);
    free(N_upper);
  }
}

#include "spamm.h"
#include "spamm_types_private.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
  for(dim = 0; dim < number_dimensions; dim++)
  {
    /* Make sure we pad at least to the extend that we can store that linear
     * kernel matrix. */
    if(use_linear_tree)
    {
      if(N[dim] < SPAMM_N_KERNEL)
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
    if(x_N > x)
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
  if(depth >= 1)
  {
    for(dim = 0; dim < number_dimensions; dim++)
    {
      if((int) (ipow(2, depth-1)) < N[dim])
      {
        depth++;
        break;
      }
    }
    depth--;
  }

  if(use_linear_tree)
  {
    while((int) (ipow(2, depth)) < SPAMM_N_KERNEL)
    {
      depth++;
    }
  }

  return depth;
}

/** Allocate a SpAMM data chunk.
 *
 * The chunk contains the following data fields. In order to guarantee this
 * layout we allocate a larger chunk of memory and then manage the data inside
 * of it ourselves.  In order to simplify access to the fields, we start the
 * chunk with a pointer array that points to the field variables.
 *
 * \code
 * struct spamm_chunk_t
 * {
 *   unsigned int *number_dimensions_pointer;
 *   unsigned int *N_block_pointer;
 *   unsigned int *N_lower_pointer;
 *   unsigned int *N_upper_pointer;
 *   float        *A_pointer;
 *   float        *A_dilated_pointer;
 *   float        *norm_pointer;
 *   float        *norm2_pointer;
 *
 *   unsigned int number_dimensions;
 *   unsigned int N_block;
 *   unsigned int N_lower[number_dimensions];
 *   unsigned int N_upper[number_dimensions];
 *
 *   spamm_float_t *A;
 *
 *   spamm_float_t *A_dilated;
 *
 *   spamm_float_t norm[];
 *   spamm_float_t norm2[];
 * };
 * \endcode
 *
 * @param number_dimensions The number of dimensions.
 * @param use_linear_tree Whether to use the linear code for the chunk or not.
 * @param N The size of original matrix (unpadded).
 * @param N_lower The lower bounds of the bounding box.
 * @param N_lower The upper bounds of the bounding box.
 *
 * @return A pointer to the newly allocated chunk.
 */
spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const short use_linear_tree,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper)
{
  void **pointer_pointer;
  unsigned int *int_pointer;

  unsigned int number_tiers;

  unsigned int *N_pointer;
  unsigned int *N_lower_pointer;
  unsigned int *N_upper_pointer;
  float *A_pointer;
  float *A_dilated_pointer;
  float *norm_pointer;
  float *norm2_pointer;

  int dim;

  spamm_chunk_t *chunk;

  chunk = spamm_allocate(spamm_chunk_get_size(number_dimensions,
        use_linear_tree, &number_tiers, N_lower, N_upper, &N_pointer,
        &N_lower_pointer, &N_upper_pointer, &A_pointer, &A_dilated_pointer,
        &norm_pointer, &norm2_pointer), 1);

  int_pointer = chunk;
  pointer_pointer = (void**) ((intptr_t) chunk + 4*sizeof(unsigned int));

  int_pointer[0] = number_dimensions;
  int_pointer[1] = number_tiers;
  int_pointer[2] = use_linear_tree;

  pointer_pointer[0] = (void*) N_pointer;
  pointer_pointer[1] = (void*) N_lower_pointer;
  pointer_pointer[2] = (void*) N_upper_pointer;
  pointer_pointer[3] = (void*) A_pointer;
  pointer_pointer[4] = (void*) A_dilated_pointer;
  pointer_pointer[5] = (void*) norm_pointer;
  pointer_pointer[6] = (void*) norm2_pointer;

  /* Store bounding box. */
  N_pointer       = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_pointer);
  N_lower_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_lower_pointer);
  N_upper_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_upper_pointer);

  for(dim = 0; dim < number_dimensions; dim++)
  {
    N_pointer[dim] = N[dim];
    N_lower_pointer[dim] = N_lower[dim];
    N_upper_pointer[dim] = N_upper[dim];
  }

  return chunk;
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

#ifdef _OPENMP
  omp_init_lock(&node->lock);
#endif

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

  for(dim = 0; dim < number_dimensions; dim++)
  {
    if(N[dim] == 0)
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
  for(dim = 0; dim < number_dimensions; dim++)
  {
    A->N[dim] = N[dim];
  }

  if(number_dimensions == 2 && use_linear_tree)
  {
    A->use_linear_tree = use_linear_tree;
  }

  /* Get tree depth. */
  A->depth = spamm_get_tree_depth(number_dimensions, A->N, use_linear_tree);

  /* Set padded matrix size. */
  A->N_padded = (int) (ipow(2, A->depth));

  /* Adjust the depth. */
  if(number_dimensions == 2 && use_linear_tree)
  {
    if(A->depth >= 4)
    {
      A->depth -= 4; /* 16x16 submatrix blocks for linear kernel. */
    }

    else
    {
      A->depth = 0;
    }
  }

  if(chunk_tier > A->depth)
  {
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

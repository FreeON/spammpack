#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Get an element from a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element \f$A(i,j)\f$.
 */
float
spamm_hashed_get (const unsigned int i, const unsigned int j, const struct spamm_hashed_t *A)
{
  unsigned int index, i_tier, j_tier;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_hashed_data_t *data;

  i_tier = (i-A->M_lower)/SPAMM_N_KERNEL;
  j_tier = (j-A->N_lower)/SPAMM_N_KERNEL;

  /* Construct linear index of the node on this tier. */
  index = spamm_index_2D(i_tier, j_tier);

  /* Get hash table at this tier. */
  node_hashtable = A->tier_hashtable[A->kernel_tier-A->tier];

  if ((data = spamm_hashtable_lookup(node_hashtable, index)) != NULL)
  {
    return data->block_dense[spamm_index_kernel_block((i-A->M_lower)%SPAMM_N_KERNEL, (j-A->N_lower)%SPAMM_N_KERNEL, A->layout)];
  }

  else { return 0; }
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
    const unsigned int contiguous_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const struct spamm_recursive_node_t *node)
{
  int dim;

  unsigned int number_rows;
  unsigned int number_columns;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  float *A;

  if (node == NULL) { return 0; }

  number_rows = N_upper[0]-N_lower[0];

  if (number_dimensions == 2 && tier == contiguous_tier && use_linear_tree)
  {
    return spamm_hashed_get(i[0], i[1], node->tree.hashed_tree);
  }

  else if (tier == contiguous_tier)
  {
    /* Get the matrix element. */
    if (node->tree.chunk == NULL) { return 0.0; }
    else
    {
      A = spamm_chunk_get_matrix(node->tree.chunk);
      return A[spamm_chunk_matrix_index(number_dimensions, N_block, N_lower, i)];
    }
  }

  else
  {
    new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < number_dimensions; dim++)
    {
      new_N_lower[dim] = N_lower[dim];
      new_N_upper[dim] = N_upper[dim];
    }

    switch (number_dimensions)
    {
      case 1:
        if (i[0] < N_lower[0]+(number_rows)/2)
        {
          new_N_upper[0] = N_lower[0]+(number_rows)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[0]);
        }

        else
        {
          new_N_lower[0] = N_lower[0]+(number_rows)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[1]);
        }

      case 2:
        number_columns = N_upper[1]-N_lower[1];
        if (i[0] < N_lower[0]+(number_rows)/2 &&
            i[1] < N_lower[1]+(number_columns)/2)
        {
          new_N_upper[0] = N_lower[0]+(number_rows)/2;
          new_N_upper[1] = N_lower[1]+(number_columns)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[0]);
        }

        else if (i[0] <  N_lower[0]+(number_rows)/2 &&
            i[1] >= N_lower[1]+(number_columns)/2)
        {
          new_N_upper[0] = N_lower[0]+(number_rows)/2;
          new_N_lower[1] = N_lower[1]+(number_columns)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[1]);
        }

        else if (i[0] >= N_lower[0]+(number_rows)/2 &&
            i[1] <  N_lower[1]+(number_columns)/2)
        {
          new_N_lower[0] = N_lower[0]+(number_rows)/2;
          new_N_upper[1] = N_lower[1]+(number_columns)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[2]);
        }

        else
        {
          new_N_lower[0] = N_lower[0]+(number_rows)/2;
          new_N_lower[1] = N_lower[1]+(number_columns)/2;
          return spamm_recursive_get(number_dimensions, i, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, N_block, use_linear_tree,
              node->tree.child[3]);
        }
        break;

      default:
        SPAMM_FATAL("not implemented\n");
    }

    free(new_N_lower);
    free(new_N_upper);
  }

  SPAMM_FATAL("I should not be here\n");
  return 0;
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
  int dim;

  unsigned int *N_lower;
  unsigned int *N_upper;

  assert(A != NULL);

  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    if (i[dim] >= A->N[dim])
    {
      SPAMM_FATAL("i[%u] out of bounds (i[%u] = %i and N[%u] = %i)\n", dim, dim, i, dim, A->N[dim]);
    }
  }

  if (A->number_dimensions == 2 && A->contiguous_tier == 0 && A->use_linear_tree)
  {
    /* In case we only have a linear tree. */
    return spamm_hashed_get(i[0], i[1], A->tree.hashed_tree);
  }

  else
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    return spamm_recursive_get(A->number_dimensions, i, N_lower, N_upper, 0,
        A->contiguous_tier, A->N_block, A->use_linear_tree,
        A->tree.recursive_tree);

    free(N_lower);
    free(N_upper);
  }
}

/** Get the number of rows of a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of rows.
 */
unsigned int
spamm_get_number_of_rows (const struct spamm_hashed_t *const A)
{
  return A->M_upper-A->M_lower;
}

/** Get the number of columns of a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of rows.
 */
unsigned int
spamm_get_number_of_columns (const struct spamm_hashed_t *const A)
{
  return A->N_upper-A->N_lower;
}

/** Get the Frobenius norm of the matrix.
 *
 * @param A The matrix.
 *
 * @return The Frobenius norm.
 */
float
spamm_hashed_get_norm (const struct spamm_hashed_t *const A)
{
  struct spamm_hashtable_t *tier_hashtable;
  struct spamm_hashed_node_t *root;

  assert(A != NULL);

  if ((tier_hashtable = A->tier_hashtable[0]) == NULL)
  {
    return 0;
  }

  if ((root = spamm_hashtable_lookup(tier_hashtable, 0)) == NULL)
  {
    return 0;
  }

  return root->norm;
}

/** Get the Frobenius norm of the matrix.
 *
 * @param A The matrix.
 *
 * @return The Frobenius norm.
 */
float
spamm_get_norm (const struct spamm_matrix_t *const A)
{
  assert(A != NULL);

  if (A->tree.recursive_tree != NULL)
  {
    return A->tree.recursive_tree->norm;
  }

  else if (A->tree.hashed_tree != NULL)
  {
    return spamm_hashed_get_norm(A->tree.hashed_tree);
  }

  else { return 0; }
}

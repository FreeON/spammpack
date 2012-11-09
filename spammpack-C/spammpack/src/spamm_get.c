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
    const unsigned int chunk_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const unsigned int depth,
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
        tier+1, chunk_tier, N_block, use_linear_tree, depth,
        node->tree.child[child_index]);

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
        A->chunk_tier, A->N_block, A->use_linear_tree, A->depth,
        A->tree.recursive_tree);

    free(N_lower);
    free(N_upper);
  }

  return Aij;
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

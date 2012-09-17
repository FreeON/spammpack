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
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element Aij.
 */
float
spamm_recursive_get (const unsigned int i, const unsigned int j, const struct spamm_recursive_node_t *node)
{
  unsigned int number_rows;
  unsigned int number_columns;

  number_rows = node->M_upper-node->M_lower;
  number_columns = node->N_upper-node->N_lower;

  if (node == NULL) { return 0; }

  else if (number_rows == node->N_linear)
  {
    return spamm_hashed_get(i, j, node->hashed_tree);
  }

  else if (number_rows == node->N_contiguous)
  {
    /* Get the matrix element. */
    if (node->data == NULL) { return 0.0; }
    else
    {
      return node->data[spamm_index_column_major(i-node->M_lower, j-node->N_lower, node->N_contiguous, node->N_contiguous)];
    }
  }

  else
  {
    if (i < node->M_lower+(number_rows)/2 &&
        j < node->N_lower+(number_columns)/2)
    {
      return spamm_recursive_get(i, j, node->child[0]);
    }

    else if (i <  node->M_lower+(number_rows)/2 &&
        j >= node->N_lower+(number_columns)/2)
    {
      return spamm_recursive_get(i, j, node->child[1]);
    }

    else if (i >= node->M_lower+(number_rows)/2 &&
        j <  node->N_lower+(number_columns)/2)
    {
      return spamm_recursive_get(i, j, node->child[2]);
    }

    else if (i >= node->M_lower+(number_rows)/2 &&
        j >= node->N_lower+(number_columns)/2)
    {
      return spamm_recursive_get(i, j, node->child[3]);
    }

    else
    {
      SPAMM_FATAL("should not be here...\n");

      /* Appease the compiler. */
      return 0;
    }
  }
}

/** Get an element from a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element.
 */
float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_matrix_t *A)
{
  assert(A != NULL);

  if (i >= A->M)
  {
    SPAMM_FATAL("i out of bounds (i = %i and M = %i)\n", i, A->M);
  }

  if (j >= A->N)
  {
    SPAMM_FATAL("j out of bounds (j = %i and N = %i)\n", j, A->N);
  }

  if (A->linear_tier == 0)
  {
    /* In case we only have a linear tree. */
    return spamm_hashed_get(i, j, A->hashed_tree);
  }

  else
  {
    return spamm_recursive_get(i, j, A->recursive_tree);
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
  return A->M;
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
  return A->N;
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

  if (A->recursive_tree != NULL)
  {
    return A->recursive_tree->norm;
  }

  else if (A->hashed_tree != NULL)
  {
    return spamm_hashed_get_norm(A->hashed_tree);
  }

  else { return 0; }
}

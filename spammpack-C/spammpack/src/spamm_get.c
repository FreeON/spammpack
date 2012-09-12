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
  unsigned int index, i_tier, j_tier, delta_index;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_hashed_data_t *data;
  float Aij = 0;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    fprintf(stderr, "illegal index values for A_ij\n");
    exit(1);
  }

  /* Go into kernel tier hash and retrieve proper node. */
  delta_index = (unsigned int) floor(A->N_padded/pow(2, A->kernel_tier));

  i_tier = i/delta_index;
  j_tier = j/delta_index;

  /* Construct linear index of the node on this tier. */
  index = spamm_index_2D(i_tier, j_tier);

  /* Get hash table at this tier. */
  node_hashtable = A->tier_hashtable[A->kernel_tier];

  if ((data = spamm_hashtable_lookup(node_hashtable, index)) != NULL)
  {
    Aij = data->block_dense[spamm_index_kernel_block(i%delta_index, j%delta_index, A->layout)];
  }

  return Aij;
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
double
spamm_recursive_get (const unsigned int i, const unsigned int j, const struct spamm_recursive_t *A)
{
  struct spamm_recursive_node_t **node = NULL;

  if (A == NULL)
  {
    return 0;
  }

  node = (struct spamm_recursive_node_t**) &(A->root);

  while (1)
  {
    if (*node == NULL)
    {
      return 0;
    }

    if ((*node)->M_upper-(*node)->M_lower == A->N_contiguous)
    {
      /* Get the matrix element. */
      if ((*node)->data == NULL) { return 0.0; }
      else
      {
        return (*node)->data[spamm_index_column_major(i-(*node)->M_lower, j-(*node)->N_lower, A->N_contiguous, A->N_contiguous)];
      }
    }

    else
    {
      if (i < (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
          j < (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[0]);
      }

      else if (i <  (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[1]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j <  (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[2]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[3]);
      }
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
spamm_get_norm (const struct spamm_hashed_t *const A)
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

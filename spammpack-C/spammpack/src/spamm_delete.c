#include "spamm.h"
#include "spamm_types_private.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Delete a hashtable entry.
 *
 * @param value The hashtable value.
 * @param user_data A pointer to an int indicating whether we are at the
 * kernel tier or not.
 */
void
spamm_delete_node_hashentry (unsigned int index, void *value, void *user_data)
{
  int *at_kernel_tier = user_data;
  struct spamm_hashed_node_t *node;
  struct spamm_hashed_data_t *data;

  if ((*at_kernel_tier) == 1)
  {
    data = value;
    spamm_delete_block(&data);
  }

  else
  {
    node = value;
    spamm_hashed_delete_node(&node);
  }
}

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_hashed_delete (struct spamm_hashed_t **A)
{
  int at_kernel_tier;
  unsigned int tier;

  /* Delete all data on all tiers. */
  for (tier = (*A)->tier; tier <= (*A)->kernel_tier; tier++)
  {
    if (tier == (*A)->kernel_tier)
    {
      at_kernel_tier = 1;
    }

    else
    {
      at_kernel_tier = 0;
    }

    /* Delete each tier hashtable. */
    spamm_hashtable_foreach((*A)->tier_hashtable[tier-(*A)->tier], spamm_delete_node_hashentry, &at_kernel_tier);

    /* Free the node hashtable. */
    spamm_hashtable_delete(&(*A)->tier_hashtable[tier-(*A)->tier]);
  }

  /* Free the tier hashtables. */
  free((*A)->tier_hashtable);

  /* Delete the matrix. */
  free(*A);
  *A = NULL;
}

/** Delete a recursive matrix.
 *
 * @param A The recursive matrix root.
 */
void
spamm_recursive_delete (struct spamm_recursive_node_t **node)
{
  int i;

  if (*node == NULL) { return; }

  for (i = 0; i < ipow(2, (*node)->number_dimensions); i++)
  {
    spamm_recursive_delete(&(*node)->child[i]);
  }
  free((*node)->child);

  if ((*node)->data != NULL)
  {
    free((*node)->data);
    (*node)->data = NULL;
  }

  if ((*node)->hashed_tree != NULL)
  {
    spamm_hashed_delete(&(*node)->hashed_tree);
  }

  free((*node)->N_lower);
  free((*node)->N_upper);

  free(*node);
  *node = NULL;
}

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_matrix_t **A)
{
  if (*A == NULL) { return; }

  if ((*A)->recursive_tree != NULL)
  {
    spamm_recursive_delete(&(*A)->recursive_tree);
  }

  else if ((*A)->hashed_tree != NULL)
  {
    spamm_hashed_delete(&(*A)->hashed_tree);
  }

  free((*A)->N);

  free(*A);
  *A = NULL;
}

/** Delete a node in a matrix.
 *
 * @param node The node in the matrix to delete.
 */
void
spamm_hashed_delete_node (struct spamm_hashed_node_t **node)
{
  free(*node);
  *node = NULL;
}

/** Delete a block of type spamm_hashed_data_t.
 *
 * @param data The block of data to free.
 */
void
spamm_delete_block (struct spamm_hashed_data_t **data)
{
  free(*data);
  *data = NULL;
}

/** Delete a SpAMM chunk.
 *
 * @param chunk The chunk.
 */
void
spamm_delete_chunk (spamm_chunk_t **chunk)
{
  free(*chunk);
  *chunk = NULL;
}

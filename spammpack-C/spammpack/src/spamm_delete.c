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
spamm_recursive_delete (const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    struct spamm_recursive_node_t **node)
{
  unsigned int i;

  if (*node == NULL) { return; }

  if (number_dimensions == 2 && use_linear_tree && tier == chunk_tier)
  {
    spamm_hashed_delete(&(*node)->tree.hashed_tree);
  }

  else if (tier == chunk_tier)
  {
    free((*node)->tree.chunk);
    (*node)->tree.chunk = NULL;
  }

  else
  {
    for (i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_delete(number_dimensions, tier+1, chunk_tier, use_linear_tree, &(*node)->tree.child[i]);
    }
    free((*node)->tree.child);
    (*node)->tree.child = NULL;
  }

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

  if ((*A)->number_dimensions == 2 && (*A)->use_linear_tree && (*A)->chunk_tier == 0)
  {
    spamm_hashed_delete(&(*A)->tree.hashed_tree);
  }

  else if ((*A)->chunk_tier == 0)
  {
    free((*A)->tree.chunk);
  }

  else
  {
    spamm_recursive_delete((*A)->number_dimensions, 0, (*A)->chunk_tier, (*A)->use_linear_tree, &(*A)->tree.recursive_tree);
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

#include "spamm.h"
#include <stdlib.h>

void
spamm_delete_node_hashentry (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int *index = key;
  int *at_kernel_tier = user_data;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  if ((*at_kernel_tier) == 1)
  {
    data = value;
    spamm_delete_block(&data);
  }

  else
  {
    node = value;
    spamm_delete_node(&node);
  }

  /* Free index key. */
  free(index);
}

void
spamm_delete_tier_hashtable (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int *tier = key;
  unsigned int *kernel_tier = user_data;
  int at_kernel_tier;
  GHashTable *node_hashtable = value;

  if ((*tier) == (*kernel_tier))
  {
    at_kernel_tier = 1;
  }

  else
  {
    at_kernel_tier = 0;
  }

  /* Delete each tier hashtable. */
  g_hash_table_foreach(node_hashtable, spamm_delete_node_hashentry,  &at_kernel_tier);

  /* Free tier key. */
  free(tier);

  /* Free the node hashtable. */
  g_hash_table_destroy(node_hashtable);
}

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_t **A)
{
  /* Delete all data on all tiers. */
  g_hash_table_foreach((*A)->tier_hashtable, spamm_delete_tier_hashtable, &(*A)->kernel_tier);

  /* Delete the tier hashtable. */
  g_hash_table_destroy((*A)->tier_hashtable);

  /* Delete the matrix. */
  free(*A);
  *A = NULL;
}

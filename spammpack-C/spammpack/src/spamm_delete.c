#include "spamm.h"
#include <stdlib.h>

void
spamm_delete_node_hashentry (void *key, void *value, void *user_data)
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

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_t **A)
{
  int at_kernel_tier;
  unsigned int tier;

  /* Delete all data on all tiers. */
  for (tier = 0; tier <= (*A)->kernel_tier; tier++)
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
    spamm_hashtable_foreach((*A)->tier_hashtable[tier], spamm_delete_node_hashentry,  &at_kernel_tier);

    /* Free the node hashtable. */
    spamm_hashtable_delete(&(*A)->tier_hashtable[tier]);
  }

  /* Free the tier hashtables. */
  free((*A)->tier_hashtable);

  /* Delete the matrix. */
  free(*A);
  *A = NULL;
}

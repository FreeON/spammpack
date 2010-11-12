#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Delete a node, and recursively all children nodes.
 *
 * This is the recursive part. Call spamm_delete() instead.
 *
 * @param node The node to delete.
 */
void
spamm_delete_node (struct spamm_node_t **node)
{
  int i, j;

  assert(*node != NULL);

  if ((*node)->block_dense != NULL)
  {
    if ((*node)->tier == (*node)->kernel_tier)
    {
      spamm_free((*node)->block_dense);
      spamm_free((*node)->block_dense_dilated);
    }
    (*node)->block_dense = NULL;
    (*node)->block_dense_dilated = NULL;
  }

  for (i = 0; i < SPAMM_N_CHILD; ++i) {
    for (j = 0; j < SPAMM_N_CHILD; ++j)
    {
      if ((*node)->child[i][j] != NULL)
      {
        spamm_delete_node(&(*node)->child[i][j]);
      }
    }
  }

  spamm_free(*node);
  *node = NULL;
}

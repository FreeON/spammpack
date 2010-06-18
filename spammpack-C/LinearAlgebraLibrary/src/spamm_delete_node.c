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
    free((*node)->block_dense);
    (*node)->block_dense = NULL;
  }

  if ((*node)->linear_quadtree != NULL)
  {
    spamm_ll_delete(NULL, &(*node)->linear_quadtree);
    spamm_mm_delete(&(*node)->linear_quadtree_memory);
  }

  if ((*node)->child != NULL)
  {
    for (i = 0; i < (*node)->M_child; ++i) {
      for (j = 0; j < (*node)->N_child; ++j)
      {
        spamm_delete_node(&(*node)->child[spamm_dense_index(i, j, (*node)->M_child, (*node)->N_child)]);
      }
    }
    free((*node)->child);
    (*node)->child = NULL;
  }

  free(*node);
  *node = NULL;
}

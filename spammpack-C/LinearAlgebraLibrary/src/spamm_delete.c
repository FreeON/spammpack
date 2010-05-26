/* @file */

#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Delete a node, and recursively all children nodes.
 *
 * This is the recursive part. Call spamm_delete() instead.
 *
 * @param node The node to delete.
 */
void
spamm_delete_node (struct spamm_node_t *node)
{
  int i, j;

  if (node->block_dense != NULL)
  {
    free(node->block_dense);
    node->block_dense = NULL;
  }

  if (node->child != NULL)
  {
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_delete_node(node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
        free(node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
        node->child[spamm_dense_index(i, j, node->M_child, node->N_child)] = NULL;
      }
    }
    free(node->child);
    node->child = NULL;
  }
}

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_t *A)
{
  assert(A != NULL);

  /* Recurse down and free() nodes. */
  //spamm_log("deleting matrix\n", __FILE__, __LINE__);
  if (A->root != NULL)
  {
    spamm_delete_node(A->root);
  }

  free(A->root);
  A->root = NULL;

  A->M = 0;
  A->N = 0;

  A->M_padded = 0;
  A->N_padded = 0;

  A->M_child = 0;
  A->N_child = 0;

  A->M_block = 0;
  A->N_block = 0;
}

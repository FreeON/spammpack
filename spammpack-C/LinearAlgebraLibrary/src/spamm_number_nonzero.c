#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Count how many nonzero nodes exist in matrix tree.
 *
 * @param node The node to count.
 *
 * @returns The number of nonzero nodes of this node and below.
 */
unsigned int
spamm_number_nonzero_node (const struct spamm_node_t *node)
{
  int i, j;
  unsigned int result = 0;

  assert(node != NULL);

  //spamm_print_node(node);
  if (node->block_dense != NULL)
  {
    //LOG("counting dense block\n");
    for (i = 0; i < node->M_block; ++i) {
      for (j = 0; j < node->N_block; ++j)
      {
        if (node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] != 0.0)
        {
          result++;
        }
      }
    }
    //LOG("counted %u nonzero elements in this block\n", result);
  }

  else if (node->child != NULL)
  {
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        result += spamm_number_nonzero_node(node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  return result;
}

/** Count how many non-zero nodes are in matrix tree.
 *
 * @param A The matrix.
 *
 * @return The number of non-zero blocks in this matrix.
 */
unsigned int
spamm_number_nonzero (const struct spamm_t *A)
{
  unsigned int result = 0;

  assert(A != NULL);

  //LOG("starting to recurse\n");
  if (A->root != NULL)
  {
    result = spamm_number_nonzero_node (A->root);
  }

  return result;
}

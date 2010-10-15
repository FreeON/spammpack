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
  if (node->tier == node->tree_depth)
  {
    LOG2_DEBUG("counting dense block\n");
    for (i = 0; i < SPAMM_M_BLOCK; ++i) {
      for (j = 0; j < SPAMM_N_BLOCK; ++j)
      {
        if (node->block_dense[spamm_dense_index(i, j, SPAMM_M_BLOCK, SPAMM_N_BLOCK)] != 0.0)
        {
          result++;
        }
      }
    }
    LOG_DEBUG("counted %u nonzero elements in this block\n", result);
  }

  for (i = 0; i < SPAMM_M_CHILD; ++i) {
    for (j = 0; j < SPAMM_N_CHILD; ++j)
    {
      if (node->child[i][j] != NULL)
      {
        result += spamm_number_nonzero_node(node->child[i][j]);
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

  LOG2_DEBUG("starting to recurse\n");
  if (A->root != NULL)
  {
    result = spamm_number_nonzero_node (A->root);
  }

  return result;
}

#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

unsigned int
spamm_number_nonzero_node (const struct spamm_node_t *node)
{
  int i, j;
  unsigned int result = 0;

  assert(node != NULL);

  //spamm_print_node(node);
  if (node->block_dense != NULL)
  {
    //spamm_log("counting dense block\n", __FILE__, __LINE__);
    for (i = 0; i < node->M_block; ++i) {
      for (j = 0; j < node->N_block; ++j)
      {
        if (node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] != 0.0)
        {
          result++;
        }
      }
    }
    //spamm_log("counted %u nonzero elements in this block\n", __FILE__, __LINE__, result);
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

unsigned int
spamm_number_nonzero (const struct spamm_t *A)
{
  unsigned int result = 0;

  assert(A != NULL);

  //spamm_log("starting to recurse\n", __FILE__, __LINE__);
  if (A->root != NULL)
  {
    result = spamm_number_nonzero_node (A->root);
  }

  return result;
}

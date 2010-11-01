#include "spamm.h"
#include <stdlib.h>

/** Initialize new node of matrix.
 *
 * @return The spamm_node_t node to initialize.
 */
struct spamm_node_t *
spamm_new_node ()
{
  unsigned int i, j;
  struct spamm_node_t *node;

  node = (struct spamm_node_t*) spamm_allocate(sizeof(struct spamm_node_t));

  node->tier = 0;
  node->tree_depth = 0;
  node->linear_tier = 0;
  node->kernel_tier = 0;

  node->M_upper = 0;
  node->M_lower = 0;
  node->N_upper = 0;
  node->N_lower = 0;

  node->M_upper_kernel_tier = 0;
  node->M_lower_kernel_tier = 0;
  node->N_upper_kernel_tier = 0;
  node->N_lower_kernel_tier = 0;

  node->index = 0;

  node->previous_i = NULL;
  node->next_i = NULL;
  node->previous_j = NULL;
  node->next_j = NULL;

  node->ordering = none;

  for (i = 0; i < SPAMM_N_CHILD; i++) {
    for (j = 0; j < SPAMM_N_CHILD; j++)
    {
      node->child[i][j] = NULL;
    }
  }

  node->linear_quadtree = NULL;
  node->linear_quadtree_default_chunksize = 1*1024*1024; /* 1 MB. */
  node->linear_quadtree_memory = NULL;
  node->block_dense = NULL;
  node->block_dense_dilated = NULL;

  node->norm = 0;
  node->norm2 = 0;

  return node;
}

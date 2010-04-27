#include "spamm.h"
#include <stdlib.h>

void
spamm_new_node (struct spamm_node_t **node)
{
  *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  (*node)->M_upper = 0;
  (*node)->M_lower = 0;
  (*node)->N_upper = 0;
  (*node)->N_lower = 0;

  (*node)->M_child = 0;
  (*node)->N_child = 0;

  (*node)->M_block = 0;
  (*node)->N_block = 0;

  (*node)->threshold = 0.0;

  (*node)->index = 0;
  (*node)->ordering = none;

  (*node)->block_loaded_in_GPU = 0;

  (*node)->child = NULL;
  (*node)->block_dense = NULL;
}

#include "spamm.h"
#include <stdlib.h>

void
spamm_new_node (struct spamm_node_t **node)
{
  int i;

  *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  (*node)->M_upper = 0;
  (*node)->M_lower = 0;
  (*node)->N_upper = 0;
  (*node)->N_lower = 0;

  (*node)->M_block = 0;
  (*node)->N_block = 0;

  for (i = 0; i < 4; ++i)
  {
    (*node)->child[i] = NULL;
  }
  (*node)->block_dense = NULL;
}

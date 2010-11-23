#include "spamm.h"
#include <stdlib.h>

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D,
    const unsigned int index_3D_ik0, const unsigned int index_3D_0kj)
{
  struct spamm_node_t *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  node->tier = tier;
  node->index_2D = index_2D;
  node->index_3D_ik0 = index_3D_ik0;
  node->index_3D_0kj = index_3D_0kj;

  return node;
}

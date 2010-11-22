#include "spamm.h"
#include <stdlib.h>

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index)
{
  struct spamm_node_t *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  node->tier = tier;
  node->index = index;

  return node;
}

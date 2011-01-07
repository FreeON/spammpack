#include "spamm.h"
#include <stdlib.h>

/** Allocate a new node of a matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param index_2D The 2D linear matrix index of this node.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D)
{
  struct spamm_node_t *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  node->tier = tier;
  node->index_2D = index_2D;

  node->norm = 0.0;
  node->norm2 = 0.0;

  return node;
}

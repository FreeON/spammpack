#include "spamm.h"
#include <stdlib.h>

/** Allocate a new node of a matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param index_2D The 2D linear matrix index of this node.
 * @param index_3D_ik0 The 3D convolution index in ik0 order.
 * @param index_3D_0kj The 3D convolution index in 0kj order.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D,
    const unsigned int index_3D_ik0, const unsigned int index_3D_0kj)
{
  struct spamm_node_t *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  node->tier = tier;
  node->index_2D = index_2D;
  node->index_3D_ik0 = index_3D_ik0;
  node->index_3D_0kj = index_3D_0kj;

  node->norm = 0.0;
  node->norm2 = 0.0;

  return node;
}

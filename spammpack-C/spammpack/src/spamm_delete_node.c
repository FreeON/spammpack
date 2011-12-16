#include "spamm.h"

#include <stdlib.h>

/** Delete a node in a matrix.
 *
 * @param node The node in the matrix to delete.
 */
void
spamm_delete_node (struct spamm_node_t **node)
{
  free(*node);
  *node = NULL;
}

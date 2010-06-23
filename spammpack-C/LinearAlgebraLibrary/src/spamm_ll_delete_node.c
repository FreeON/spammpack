#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Delete an element of a linked list.
 *
 * @param node The node to delete.
 * @param list The list from which to delete the node.
 *
 * \bug This function currently does not check whether the node is in the list
 *      or not.
 */
void
spamm_ll_delete_node (struct spamm_ll_node_t *node, struct spamm_ll_t *list)
{
  assert(node != NULL);
  assert(list != NULL);

  /* Check whether the node is in the list. */

  /* Delete the node. */
  if (list->first == node)
  {
  }
}

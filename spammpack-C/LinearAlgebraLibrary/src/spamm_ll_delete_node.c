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
  if (list->number_elements == 0)
  {
    LOG2_FATAL("list has no elements\n");
    exit(1);
  }

  /* Delete the node. */
  if (list->first == node)
  {
    list->first = node->next;
  }

  if (list->last == node)
  {
    list->last = node->previous;
  }

  if (node->previous != NULL)
  {
    node->previous->next = node->next;
  }

  if (node->next != NULL)
  {
    node->next->previous = node->previous;
  }

  list->number_elements--;
}

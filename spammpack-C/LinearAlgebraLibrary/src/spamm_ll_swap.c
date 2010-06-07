#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Swap to elements in a linked list.
 *
 * @param node1 The first node of the list.
 * @param node2 The second node of the list.
 * @param list The linked list.
 */
void
spamm_ll_swap (struct spamm_ll_node_t **node1, struct spamm_ll_node_t **node2,
    struct spamm_ll_t *list)
{
  struct spamm_ll_node_t *node1_previous = (*node1)->previous;
  struct spamm_ll_node_t *node2_previous = (*node2)->previous;
  struct spamm_ll_node_t *node1_next = (*node1)->next;
  struct spamm_ll_node_t *node2_next = (*node2)->next;

  struct spamm_ll_node_t *temp;

  assert(list != NULL);
  assert(*node1 != NULL);
  assert(*node2 != NULL);

  /* Swap out first and last link. */
  if      (list->first == *node1) { list->first = *node2; }
  else if (list->first == *node2) { list->first = *node1; }
  if      (list->last == *node1)  { list->last = *node2; }
  else if (list->last == *node2)  { list->last = *node1; }

  /* Connect neighbors of node1 and node2. */
  if (node1_previous != NULL) { node1_previous->next = *node2; }
  if (node2_previous != NULL) { node2_previous->next = *node1; }
  if (node1_next != NULL) { node1_next->previous = *node2; }
  if (node2_next != NULL) { node2_next->previous = *node1; }

  /* Connect to neighbors. */
  temp = (*node1)->previous; (*node1)->previous = (*node2)->previous; (*node2)->previous = temp;
  temp = (*node1)->next; (*node1)->next = (*node2)->next; (*node2)->next = temp;

  /* Swap nodes. */
  temp = *node1; *node1 = *node2; *node2 = temp;
}

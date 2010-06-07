#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Append a node to a linked list.
 *
 * @param data The data to append.
 * @param list The linked list.
 */
void
spamm_ll_append (void *data, struct spamm_ll_t *list)
{
  struct spamm_ll_node_t *new_node = NULL;

  spamm_ll_initialize_node(&new_node);

  /* Append new node to list. */
  new_node->previous = list->last;
  if (list->first == NULL) { list->first = new_node; }
  if (list->last != NULL) { list->last->next = new_node; }
  list->last = new_node;

  /* Store data. */
  new_node->data = data;

  /* Increment element counter. */
  list->number_elements++;
}

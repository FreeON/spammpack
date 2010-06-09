#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Delete a linked list.
 *
 * @param list The linked list to delete.
 */
void
spamm_ll_delete (struct spamm_ll_t *list)
{
  struct spamm_ll_node_t *node1, *node2;

  assert(list != NULL);

  for (node1 = list->first; node1 != NULL; )
  {
    node2 = node1;
    node1 = node1->next;

    free(node2);
  }
  spamm_ll_initialize(list);
}

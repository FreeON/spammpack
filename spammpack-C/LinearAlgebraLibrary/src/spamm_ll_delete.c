#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Delete a linked list.
 *
 * @param delete_data A function that takes care of deleting the data. In a
 *                    simple case this could be free().
 * @param list The linked list to delete.
 */
void
spamm_ll_delete (void (*delete_data) (void*), struct spamm_ll_t **list)
{
  struct spamm_ll_node_t *node1, *node2;

  assert(*list != NULL);

  for (node1 = (*list)->first; node1 != NULL; )
  {
    if (delete_data != NULL)
    {
      delete_data(node1->data);
    }

    node2 = node1;
    node1 = node1->next;

    free(node2);
  }

  free(*list);
  *list = NULL;
}

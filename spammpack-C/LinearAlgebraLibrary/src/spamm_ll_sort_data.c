#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Sort a linked list.
 *
 * This sort copies the data in the linked list and does not touch the links
 * between the list elements. Compare with spamm_ll_sort() which does the
 * opposite.
 *
 * @param compare A function that compares data1 and data2 and returns -1 for
 *                data1 < data2, 0 for data1 == data2, and +1 for data1 >
 *                data2.
 * @param swap A function that swaps data1 and data2.
 * @param list The linked list.
 */
void
spamm_ll_sort_data (int (*compare) (const void *data1, const void *data2),
    void (*swap) (void *data1, void *data2), struct spamm_ll_t *list)
{
  struct spamm_ll_node_t *node1;
  struct spamm_ll_node_t *node2;

  assert(list != NULL);

  for (node1 = list->first; node1 != list->last; node1 = node1->next) {
    for (node2 = node1->next; node2 != NULL; node2 = node2->next)
    {
      if (compare(node1->data, node2->data) > 0)
      {
        swap(node1->data, node2->data);
      }
    }
  }
}

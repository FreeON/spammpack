#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Sort a linked list.
 *
 * @param compare A function that compares data1 and data2 and returns -1 for
 *                data1 < data2, 0 for data1 == data2, and +1 for data1 >
 *                data2.
 * @param list The linked list.
 */
void
spamm_ll_sort (int (*compare) (const void *data1, const void *data2), struct spamm_ll_t *list)
{
  /* We optimize locality by sorting the stream so as to minimize the change
   * of any of the indices. */
  struct spamm_ll_node_t *node1, *node2;

  assert(list != NULL);

  for (node1 = list->first; node1 != list->last && node1 != NULL; node1 = node1->next) {
    for (node2 = node1->next; node2 != NULL; node2 = node2->next)
    {
      if (compare(node1->data, node2->data) > 0)
      {
        spamm_ll_swap(&node1, &node2, list);
      }
    }
  }
}

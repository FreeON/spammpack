#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Get an element from a linked list.
 *
 * @param i Index of element to get. Counting starts with 0.
 * @param list The linked list.
 *
 * @return The element of the linked list. A return value of NULL means that
 *         the element was not found in the list.
 */
struct spamm_ll_node_t *
spamm_ll_get (const unsigned int i, const struct spamm_ll_t *list)
{
  int j;
  struct spamm_ll_node_t *node;

  assert(list != NULL);

  if (i >= list->number_elements)
  {
    return NULL;
  }

  j = 0;
  for (node = list->first; node != NULL; node = node->next)
  {
    if (i == j) { break; }
    j++;
  }

  if (j < list->number_elements && node == NULL)
  {
    LOG("bug? i = %i, j = %i\n", i, j);
    spamm_ll_print(list);
    exit(1);
  }

  return node;
}

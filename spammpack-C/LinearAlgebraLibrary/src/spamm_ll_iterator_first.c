#include "spamm.h"
#include <stdlib.h>

/** List iterator.
 *
 * @param iterator The iterator.
 *
 * @return First element of list.
 */
struct spamm_ll_node_t *
spamm_ll_iterator_first (struct spamm_ll_iterator_t *iterator)
{
  iterator->last_node = iterator->list->first;
  return iterator->last_node;
}

#include "spamm.h"
#include <stdlib.h>

/** Iterator: next.
 *
 * The iterator must be initialized with a call to spamm_ll_iterator_first().
 * Otherwise the internal state of the iterator is undefined.
 *
 * @param iterator The iterator.
 *
 * @return The next element in the list this iterator operates on.
 */
struct spamm_ll_node_t *
spamm_ll_iterator_next (struct spamm_ll_iterator_t *iterator)
{
  if (iterator->last_node != NULL)
  {
    iterator->last_node = iterator->last_node->next;
  }

  return iterator->last_node;
}

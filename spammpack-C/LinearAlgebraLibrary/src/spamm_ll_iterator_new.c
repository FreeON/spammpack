#include "spamm.h"
#include <stdlib.h>

/** Create a new iterator on a linked list.
 *
 * @param list The list to apply the iterator to.
 *
 * @return The iterator.
 */
struct spamm_ll_iterator_t *
spamm_ll_iterator_new (struct spamm_ll_t *list)
{
  struct spamm_ll_iterator_t *result;

  result = (struct spamm_ll_iterator_t*) malloc(sizeof(struct spamm_ll_iterator_t));

  if (result != NULL)
  {
    result->list = list;
    result->last_node = NULL;
  }

  else
  {
    spamm_log("failed to allocate iterator\n", __FILE__, __LINE__);
  }

  return result;
}

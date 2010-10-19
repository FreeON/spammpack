#include "spamm.h"
#include <stdlib.h>

/** Delete an iterator.
 *
 * @param iterator The iterator to delete.
 */
void
spamm_ll_iterator_delete (struct spamm_ll_iterator_t **iterator)
{
  free(*iterator);
  *iterator = NULL;
}

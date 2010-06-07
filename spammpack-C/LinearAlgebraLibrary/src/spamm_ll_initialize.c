#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Initialize a new linked list.
 *
 * @param list The linked list to initialize.
 */
void
spamm_ll_initialize (struct spamm_ll_t *list)
{
  assert(list != NULL);

  list->number_elements = 0;
  list->first = NULL;
  list->last = NULL;
}

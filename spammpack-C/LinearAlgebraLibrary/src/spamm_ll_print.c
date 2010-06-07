#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Print a linked list.
 *
 * @param list The linked list to print.
 */
void
spamm_ll_print (const struct spamm_ll_t *list)
{
  unsigned int i;

  assert(list != NULL);

  struct spamm_ll_node_t *node;

  printf("linked list: [ ");
  i = 0;
  for (node = list->first; node != NULL; node = node->next)
  {
    printf("%u:%p ", ++i, node->data);
  }
  printf("]\n");
}

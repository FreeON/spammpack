#include "spamm.h"
#include "spamm_ll.h"
#include <assert.h>
#include <stdlib.h>

/** Print a linked list.
 *
 * @param data_to_string A function that converts data into a string.
 * @param list The linked list to print.
 */
void
spamm_ll_print (char *(*data_to_string) (const void *data), const struct spamm_ll_t *list)
{
  char *data_string;
  unsigned int i;

  assert(list != NULL);

  struct spamm_ll_node_t *node;

  printf("linked list: [ ");
  i = 0;
  for (node = list->first; node != NULL; node = node->next)
  {
    data_string = data_to_string(node->data);
    printf("%u:%p (%s) ", ++i, node->data, data_string);
    free(data_string);
  }
  printf("]\n");
}

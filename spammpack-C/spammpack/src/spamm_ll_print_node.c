#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Print information of a linked list node.
 *
 * @param data_to_string A function that converts data into a string.
 * @param node The linked list node.
 */
void
spamm_ll_print_node (char *(*data_to_string) (const void *data), const struct spamm_ll_node_t *node)
{
  char *datastring;

  assert(node != NULL);

  datastring = data_to_string(node->data);
  printf("%s ", datastring);
  free(datastring);
}

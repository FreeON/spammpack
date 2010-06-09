#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Initialize a new node of a linked list.
 *
 * @param node The linked list node to initialize.
 */
void
spamm_ll_initialize_node (struct spamm_ll_node_t **node)
{
  assert(node != NULL);
  assert(*node == NULL);

  *node = (struct spamm_ll_node_t*) malloc(sizeof(struct spamm_ll_node_t));

  (*node)->previous = NULL;
  (*node)->next = NULL;
  (*node)->data = NULL;
}

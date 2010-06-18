#include "spamm.h"
#include <stdlib.h>

/** Initialize a new node of a linked list.
 *
 * @return A pointer to the newly created and initialized node.
 */
struct spamm_ll_node_t *
spamm_ll_new_node ()
{
  struct spamm_ll_node_t *result;

  result = (struct spamm_ll_node_t*) malloc(sizeof(struct spamm_ll_node_t));

  if (result != NULL)
  {
    result->previous = NULL;
    result->next = NULL;
    result->data = NULL;
  }

  else
  {
    LOG2("failed to allocated new node\n");
  }

  return result;
}

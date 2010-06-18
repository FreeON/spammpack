#include "spamm.h"
#include <stdlib.h>

/** Initialize a new linked list.
 *
 * This function allocates a new linked list and initializes it. The caller
 * has to take care of free()'ing this list at some point by calling
 * spamm_ll_delete().
 *
 * @return A pointer to the newly allocated list.
 */
struct spamm_ll_t *
spamm_ll_new ()
{
  struct spamm_ll_t *result;

  result = (struct spamm_ll_t*) malloc(sizeof(struct spamm_ll_t));

  if (result != NULL)
  {
    result->number_elements = 0;
    result->first = NULL;
    result->last = NULL;
  }

  else
  {
    LOG2_FATAL("error allocating new list\n");
  }

  return result;
}

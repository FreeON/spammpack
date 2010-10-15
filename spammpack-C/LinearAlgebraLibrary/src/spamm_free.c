#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Free a previously allocated memory chunk.
 */
void
spamm_free (void *data)
{
  assert(data != NULL);

  LOG_DEBUG("freeing memory at %p\n", data);
  free(data);
}

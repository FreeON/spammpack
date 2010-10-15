#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Free a previously allocated memory chunk.
 *
 * @param data The pointer to the block to be free()'ed.
 */
void
spamm_free (void *data)
{
  assert(data != NULL);

  LOG_DEBUG("freeing memory at %p\n", data);
  free(data);
}

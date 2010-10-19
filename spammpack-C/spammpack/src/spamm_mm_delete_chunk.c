#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Delete a chunk.
 *
 * @param data The chunk to delete.
 */
void
spamm_mm_delete_chunk (void *data)
{
  struct spamm_mm_chunk_t *chunk = data;

  assert(data != NULL);

  if (chunk->allocated_start != NULL)
  {
    spamm_ll_delete(NULL, &chunk->allocated_start);
  }

  if (chunk->allocated_end != NULL)
  {
    spamm_ll_delete(NULL, &chunk->allocated_end);
  }

  if (chunk->data != NULL)
  {
    free(chunk->data);
  }

  free(data);
}

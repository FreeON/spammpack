#include "spamm.h"
#include <stdlib.h>

/** Grow memory by one chunk.
 *
 * @param chunksize The size of contiguous chunks.
 * @param memory The linked list of memory chunks.
 */
void *
spamm_mm_grow (const unsigned int chunksize, struct spamm_mm_t *memory)
{
  struct spamm_mm_chunk_t *chunk;

  chunk = spamm_mm_new_chunk(chunksize);

  if (chunk != NULL)
  {
    spamm_ll_append(chunk, memory->chunks);
#ifdef SPAMM_MM_DEBUG
    LOG("allocated chunk of %u bytes from %p to %p\n", chunksize, chunk->data, ((char*) chunk->data)+chunksize);
#endif
  }

  else
  {
    spamm_log("failed to allocate chunk\n", __FILE__, __LINE__);
  }

  return chunk;
}

#include "spamm.h"
#include <stdlib.h>

/** Grow memory by one chunk.
 *
 * @param chunksize The size of contiguous chunks.
 * @param memory The linked list of memory chunks.
 */
void *
spamm_mm_grow (const unsigned int chunksize, struct spamm_ll_t *memory)
{
  struct spamm_mm_t *chunk = NULL;

  chunk = (struct spamm_mm_t*) malloc(sizeof(struct spamm_mm_t));

  if (chunk != NULL)
  {
    spamm_ll_append(chunk, memory);
    chunk->chunksize = chunksize;
    chunk->bytes_allocated = 0;

    /* Reset list of allocated pointers. */
    spamm_ll_initialize(&chunk->allocated_start);
    spamm_ll_initialize(&chunk->allocated_end);

    /* Allocate memory for data in chunk. */
    chunk->data = malloc(chunksize);
    if (chunk->data == NULL)
    {
      spamm_log("failed to allocate data\n", __FILE__, __LINE__);
      return NULL;
    }

    //LOG("allocated chunk of %u bytes at %p --> %p\n", chunksize, chunk->data, ((char*) chunk->data)+chunksize);
  }

  else
  {
    spamm_log("failed to allocate chunk\n", __FILE__, __LINE__);
    return NULL;
  }

  return chunk;
}

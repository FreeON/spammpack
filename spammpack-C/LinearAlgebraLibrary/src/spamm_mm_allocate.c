#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Allocate memory in spamm_mm_t managed memory.
 *
 * @param size The size of the requested memory in bytes.
 * @param memory The allocated spamm_mm_t memory.
 */
void *
spamm_mm_allocate (const unsigned int size, struct spamm_mm_t *memory)
{
  unsigned int chunksize;
  void *result = NULL;
  struct spamm_mm_chunk_t *chunk;

  assert(memory != NULL);

#ifdef SPAMM_MM_DEBUG
  spamm_log("starting...\n", __FILE__, __LINE__);
  spamm_mm_print(memory);
#endif

  chunksize = memory->chunksize;
#ifdef SPAMM_MM_DEBUG
  LOG("memory has a chunksize of %u bytes\n", chunksize);
#endif

  if (size > chunksize)
  {
    spamm_log("requested size is larger than chunksize\n", __FILE__, __LINE__);
  }

  else
  {
    chunk = ((struct spamm_mm_chunk_t*) memory->chunks->last->data);
#ifdef SPAMM_MM_DEBUG
    LOG("last memory chunk at %p with data at %p\n", chunk, chunk->data);
#endif

    if (chunk->bytes_allocated+size > chunksize)
    {
      /* Data will not fit. Allocate new chunk. */
#ifdef SPAMM_MM_DEBUG
      LOG("not enough room to fit %u bytes in last chunk. Allocating new chunk\n", size);
#endif
      spamm_mm_grow(chunksize, memory);
      chunk = ((struct spamm_mm_chunk_t*) memory->chunks->last->data);
    }

    /* Data will fit into this chunk. */
#ifdef SPAMM_MM_DEBUG
    LOG("current chunk allocation: %u bytes (out of %u bytes)\n", chunk->bytes_allocated, chunksize);
#endif

    /* Mark region in chunk as allocated. */
    if (chunk->allocated_start->number_elements == 0)
    {
#ifdef SPAMM_MM_DEBUG
      spamm_log("allocating first memory block\n", __FILE__, __LINE__);
#endif
      spamm_ll_append(((char*) chunk->data), chunk->allocated_start);
      spamm_ll_append(((char*) chunk->data)+size-1, chunk->allocated_end);
    }

    else
    {
#ifdef SPAMM_MM_DEBUG
      spamm_log("allocating memory block\n", __FILE__, __LINE__);
#endif
      spamm_ll_append(((char*) chunk->allocated_end->last->data)+1, chunk->allocated_start);
      spamm_ll_append(((char*) chunk->allocated_end->last->data)+1+size-1, chunk->allocated_end);
    }

#ifdef SPAMM_MM_DEBUG
    LOG("allocating region of %u bytes from %p to %p\n", size, chunk->allocated_start->last->data, chunk->allocated_end->last->data);
#endif
    chunk->bytes_allocated += size;
    result = chunk->allocated_start->last->data;
  }

#ifdef SPAMM_MM_DEBUG
  spamm_log("done.\n", __FILE__, __LINE__);
  spamm_mm_print(memory);
#endif

  /* Return pointer to allocated memory. */
  return result;
}

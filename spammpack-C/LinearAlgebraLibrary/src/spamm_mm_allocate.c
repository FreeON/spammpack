#include "spamm.h"
#include "spamm_ll.h"
#include "spamm_mm.h"
#include <assert.h>
#include <stdlib.h>

/** Allocate a chunk of memory.
 *
 * @param size The size of the requested memory in bytes.
 * @param memory The allocated spamm_mm_t memory.
 */
void *
spamm_mm_allocate (const unsigned int size, struct spamm_ll_t *memory)
{
  unsigned int chunksize;
  void *result = NULL;
  struct spamm_mm_t *chunk;

  assert(memory != NULL);

  chunksize = ((struct spamm_mm_t*) memory->first->data)->chunksize;
  //LOG("memory has a chunksize of %u bytes\n", chunksize);

  if (size > chunksize)
  {
    spamm_log("requested size is larger than chunksize\n", __FILE__, __LINE__);
  }

  else
  {
    chunk = ((struct spamm_mm_t*) memory->last->data);

    if (chunk->bytes_allocated+size > chunksize)
    {
      /* Data will not fit. Allocate new chunk. */
      //LOG("not enough room to fit %u bytes in last chunk. Allocating new chunk\n", size);
      spamm_mm_grow(chunksize, memory);
      chunk = ((struct spamm_mm_t*) memory->last->data);
    }

    /* Data will fit into this chunk. */
    //LOG("current chunk allocation: %u bytes (out of %u bytes)\n", ((struct spamm_mm_t*) memory->last->data)->bytes_allocated, chunksize);

    /* Mark region in chunk as allocated. */
    if (chunk->allocated_start.number_elements == 0)
    {
      //spamm_log("allocating first chunk\n", __FILE__, __LINE__);
      spamm_ll_append(((char*) chunk->data), &chunk->allocated_start);
      spamm_ll_append(((char*) chunk->data)+size-1, &chunk->allocated_end);
    }

    else
    {
      spamm_ll_append(((char*) chunk->allocated_end.last->data)+1, &chunk->allocated_start);
      spamm_ll_append(((char*) chunk->allocated_end.last->data)+1+size-1, &chunk->allocated_end);
    }

    //LOG("allocating region of %u bytes from %p to %p\n", size, chunk->allocated_start.last->data, chunk->allocated_end.last->data);
    chunk->bytes_allocated += size;
    result = chunk->allocated_start.last->data;
  }

  /* Return pointer to allocated memory. */
  //spamm_mm_print(memory);
  return result;
}

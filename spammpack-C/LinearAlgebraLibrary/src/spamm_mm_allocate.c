#include "spamm.h"
#include "spamm_mm.h"
#include <stdlib.h>

/** Allocate a chunk of memory.
 *
 * @param size The size of the requested memory in bytes.
 * @param memory The allocated spamm_mm_t memory.
 */
void *
spamm_mm_allocate (const unsigned int size, struct spamm_mm_t *memory)
{
  void *result = NULL;

  if (size > memory->chunksize)
  {
    spamm_log("[FIXME] requested memory size is larger than chunksize\n", __FILE__, __LINE__);
    exit(1);
  }

  if (memory->bytes_allocated+size <= memory->chunksize)
  {
    /* Data will fit into this chunk. */
    LOG("current chunk allocation: %u bytes (out of %u bytes)\n", memory->bytes_allocated, memory->chunksize);
  }

  else
  {
    /* Data will not fit. Allocate new chunk. */
  }

  /* [FIXME] */
  result = malloc(size);

  /* Return pointer to allocated memory. */
  return result;
}

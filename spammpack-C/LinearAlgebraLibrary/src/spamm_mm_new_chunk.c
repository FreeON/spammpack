#include "spamm.h"
#include <stdlib.h>

/** Allocate a new memory chunk.
 *
 * @param chunksize The chunksize in bytes.
 *
 * @return The new memory chunk.
 */
struct spamm_mm_chunk_t *
spamm_mm_new_chunk (const unsigned int chunksize)
{
  struct spamm_mm_chunk_t *result;

  result = (struct spamm_mm_chunk_t*) malloc(sizeof(struct spamm_mm_chunk_t));

  if (result != NULL)
  {
    result->chunksize = chunksize;
    result->bytes_allocated = 0;
    result->allocated_start = spamm_ll_new();
    result->allocated_end = spamm_ll_new();
    result->data = malloc(chunksize);
  }

  else
  {
    spamm_log("error allocating new chunk\n", __FILE__, __LINE__);
  }

  return result;
}

#include "spamm.h"
#include "spamm_mm.h"
#include <stdlib.h>

/** Initialize a new dynamic memory region.
 *
 * @param chunksize The size of contigueous chunks.
 *
 * @return A pointer to the allocated memory. A return value of NULL indicates
 *         that allocation failed.
 */
struct spamm_mm_t *
spamm_mm_initialize (const unsigned int chunksize)
{
  struct spamm_mm_t *result;

  LOG("allocating chunk of %u bytes\n", chunksize);
  result = (struct spamm_mm_t*) malloc(sizeof(struct spamm_mm_t));

  if (result != NULL)
  {
    result->next = NULL;
    result->previous = NULL;
    result->chunksize = chunksize;
    result->bytes_allocated = 0;
    result->data = NULL;
  }

  return result;
}

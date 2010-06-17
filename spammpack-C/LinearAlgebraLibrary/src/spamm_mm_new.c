#include "spamm.h"
#include <stdlib.h>

/** Initialize a new dynamic memory region.
 *
 * The function allocates spamm_mm_t and returns a pointer to it. The caller
 * needs to take care of free()'ing this memory at some point by calling
 * spamm_mm_delete().
 *
 * @param chunksize The size of contiguous chunks.
 *
 * @return A pointer to a linked list that contains the allocated memory. A
 *         return value of NULL indicates that allocation failed.
 */
struct spamm_mm_t *
spamm_mm_new (const unsigned int chunksize)
{
  struct spamm_mm_t *result;

#ifdef SPAMM_MM_DEBUG
  LOG("allocating memory with chunks of %u bytes\n", chunksize);
#endif
  result = (struct spamm_mm_t*) malloc(sizeof(struct spamm_mm_t));
  if (result != NULL)
  {
    result->chunksize = chunksize;
    result->chunks = spamm_ll_new();
    if (result->chunks != NULL)
    {
      /* Allocate memory for data chunk bookkeeping structure. */
      if (spamm_mm_grow(chunksize, result) == NULL)
      {
        spamm_log("failed to grow memory\n", __FILE__, __LINE__);
        return NULL;
      }
    }

    else
    {
      /* Allocation error occurred while trying to get memory for the linked
       * list. */
      spamm_log("allocation of linked list failed\n", __FILE__, __LINE__);
      free(result);
      result = NULL;
    }
  }

  /* Debug. */
#ifdef SPAMM_MM_DEBUG
  spamm_mm_print(result);
#endif

  return result;
}

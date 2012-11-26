#include "config.h"
#include "spamm.h"

#ifdef HAVE_POSIX_MEMALIGN
#include <errno.h>
#include <stdio.h>
#endif

#include <stdlib.h>
#include <string.h>

/** Allocate a contiguous memory block.
 *
 * The function spamm_allocate() allocates a contiguous memory block. It
 * basically uses malloc() to achieve this and tries to align the allocated
 * block on the given boundary (controlled by SPAMM_ALIGNMENT).
 *
 * @param size The size in bytes of the requested memory block.
 * @param zero_memory If set to zero, then nothing is done to the newly
 * allocated memory. The contents are undefined. If set to something else,
 * then the newly allocate memory is set to zero.
 *
 * @return A pointer to the allocated memory. A NULL pointer indicates that
 * something went wrong.
 */
void *
spamm_allocate (const size_t size, const short zero_memory)
{
  void *data;

#ifdef HAVE_POSIX_MEMALIGN
  int result;

  if((result = posix_memalign(&data, SPAMM_ALIGNMENT, size)) != 0)
  {
    switch(result)
    {
      case EINVAL:
        SPAMM_FATAL("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        SPAMM_FATAL("There was insufficient memory to fulfill the allocation request.\n");
        break;

      default:
        SPAMM_FATAL("unknown error code\n");
        break;
    }
  }

  else
  {
    /* Check whether we should zero the new memory. */
    if(zero_memory != 0)
    {
      memset(data, 0, size);
    }

    return data;
  }
#else
  /* In order to allocate aligned memory, we could allocate a larger chunk and
   * then calculate a pointer into the chunk which is correctly aligned. This
   * would necessitate our own memory management since the returned, aligned
   * pointer is not useful anymore for a free().
   */
  SPAMM_WARN("FIXME: can not allocate aligned memory\n");
  data = calloc(1, size);
  return data;
#endif

  SPAMM_FATAL("I should not be here\n");
  return NULL;
}

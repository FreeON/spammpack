#include "config.h"
#include "spamm.h"

#ifdef HAVE_POSIX_MEMALIGN
#include <errno.h>
#include <stdio.h>
#endif

#include <stdlib.h>

/** Allocate a contiguous memory block.
 *
 * The function spamm_allocate() allocates a contiguous memory block. It
 * basically uses malloc() to achieve this and tries to align the allocated
 * block on the given boundary (controlled by SPAMM_ALIGNMENT).
 *
 * @param size The size in bytes of the requested memory block.
 *
 * @return A pointer to the allocated memory. A NULL pointer indicates that
 * something went wrong.
 */
void *
spamm_allocate (size_t size)
{
  void *data;

#ifdef HAVE_POSIX_MEMALIGN
  int result;

  if ((result = posix_memalign(&data, SPAMM_ALIGNMENT, size)) != 0)
  {
    switch (result)
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
    return data;
  }
#else
  /* In order to allocate aligned memory, we could allocate a larger chunk and
   * then calculate a pointer into the chunk which is correctly aligned. This
   * would necessitate our own memory management since the returned, aligned
   * pointer is not useful anymore for a free().
   */
  printf("[%s:%i] FIXME: can not allocate aligned memory\n", __FILE__, __LINE__);
  data = calloc(size, 1);
  return data;
#endif
}


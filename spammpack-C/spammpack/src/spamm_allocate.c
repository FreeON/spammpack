#include "spamm_config.h"
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
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        exit(1);
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        exit(1);
        break;

      default:
        printf("unknown error code\n");
        exit(1);
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
  printf("[FIXME] can not allocate aligned memory\n");
  data = malloc(size);
  return data;
#endif
}


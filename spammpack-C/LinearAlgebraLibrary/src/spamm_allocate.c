#include "config.h"
#include "spamm.h"

#include <errno.h>
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
        LOG2_FATAL("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        exit(1);
        break;

      case ENOMEM:
        LOG2_FATAL("There was insufficient memory to fulfill the allocation request.\n");
        exit(1);
        break;

      default:
        LOG2_FATAL("unknown error code\n");
        exit(1);
        break;
    }
  }

  else
  {
    LOG_DEBUG("allocated %u bytes at %p\n", size, (void*) data);
    return data;
  }
#else
  data = malloc(size);
  LOG_DEBUG("allocated %u bytes at %p\n", size, (void*) data);
  return data;
#endif
}

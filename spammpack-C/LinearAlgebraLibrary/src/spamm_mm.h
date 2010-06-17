/** @file */

#if ! defined(__SPAMM_MM_H)

/** Define in case spamm_mm.h has been included. */
#define __SPAMM_MM_H 1

/** \page page_mm The SpAMM memory manager.
 *
 * \brief A memory manager to manage contiguous memory allocation.
 *
 * In order to allow for dynamic memory allocation in a contiguous chunk of
 * memory, the SpAMM memory manager provides some basic functions. Details can
 * be found in spamm_mm.h.
 *
 * \section Introduction
 *
 * The SpAMM memory manager will allocate dynamic memory in chunks of a
 * certain size. At first, only one chunk is allocated. As more memory is
 * requested, the memory manager might allocate more chunks, but links them
 * with the first one. The management of what chunk holds the data is done by
 * the memory manager.
 */

/** A memory chunk.
 */
struct spamm_mm_t
{
  /** How large is this chunk? */
  unsigned int chunksize;

  /** How many bytes are already allocated in this chunk? */
  unsigned int bytes_allocated;

  /** A list of pointers to start of "allocated" pieces of memory in the data
   * chunk. */
  struct spamm_ll_t allocated_start;

  /** A list of pointers to end of "allocated" pieces of memory in the data
   * chunk. */
  struct spamm_ll_t allocated_end;

  /** The data. */
  void *data;
};

void *
spamm_mm_allocate (const unsigned int size, struct spamm_ll_t *memory);

void *
spamm_mm_grow (const unsigned int chunksize, struct spamm_ll_t *memory);

struct spamm_ll_t *
spamm_mm_initialize (const unsigned int chunksize);

void
spamm_mm_print (const struct spamm_ll_t *memory);

#endif

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
 * \section spamm_mm_sec_introduction Introduction
 *
 * The SpAMM memory manager will allocate dynamic memory in chunks of a
 * certain size. The spamm_mm_allocate() API is just like malloc(). The
 * function returns a pointer to the start of a memory region of the requested
 * size. Different to malloc() however, spamm_mm_alloc() allocates memory in a
 * larger contiguous chunk. This ensures that several spamm_mm_allocate()
 * request are contiguous with respect to each other, resulting in packed
 * memory of chunks of chunksize bytes.
 *
 * The memory manager internally does this by first allocating a chunk of
 * memory of chunksize bytes. As this chunk is filled by calls to
 * spamm_mm_allocate(), the memory manager might allocate more chunks, but
 * links them with the first one. The management of what chunk holds the data
 * is done by the memory manager.
 *
 * \section spamm_mm_sec_typical_usage Typical usage
 *
 * First, a new memory manager object has to be created with spamm_mm_new().
 * The chunksize has to be chosen. Currently the chunksize has to be larger
 * than the largest memory block one anticipates to ever allocate. The reason
 * is simply that the chunks themselves are not allocated contiguously, and
 * hence, a block that larger than the chunksize would not fit into a
 * contiguous chunk. The whole point of the memory manager however is to
 * guarantee contiguous allocation of several memory block. Blocks are
 * allocated with a call to spamm_mm_allocate(). A pointer is returned to the
 * start of the allocated memory block.
 *
 * A typical code fragment:
 *
 * \code
 * struct spamm_mm_t *memory;
 * int *i;
 *
 * memory = spamm_mm_new(chunksize);
 * i = spamm_mm_allocate(sizeof(int), memory);
 *
 * *i = 1;
 *
 * spamm_mm_delete(&memory);
 * \endcode
 */

/** Memory managed by this memory manager. */
struct spamm_mm_t
{
  /** The chunksize of this memory. */
  unsigned int chunksize;

  /** A linked list of memory chunks. */
  struct spamm_ll_t *chunks;
};

/** A memory chunk.
 */
struct spamm_mm_chunk_t
{
  /** The chunksize of this memory. */
  unsigned int chunksize;

  /** How many bytes are already allocated in this chunk? */
  unsigned int bytes_allocated;

  /** A list of pointers to start of "allocated" pieces of memory in the data
   * chunk. */
  struct spamm_ll_t *allocated_start;

  /** A list of pointers to end of "allocated" pieces of memory in the data
   * chunk. */
  struct spamm_ll_t *allocated_end;

  /** The data. */
  void *data;
};

void *
spamm_mm_allocate (const unsigned int size, struct spamm_mm_t *memory);

void
spamm_mm_delete (struct spamm_mm_t **memory);

void
spamm_mm_delete_chunk (void *data);

void *
spamm_mm_grow (const unsigned int chunksize, struct spamm_mm_t *memory);

struct spamm_mm_t *
spamm_mm_new (const unsigned int chunksize);

struct spamm_mm_chunk_t *
spamm_mm_new_chunk (const unsigned int chunksize);

void
spamm_mm_print (const struct spamm_mm_t *memory);

#endif

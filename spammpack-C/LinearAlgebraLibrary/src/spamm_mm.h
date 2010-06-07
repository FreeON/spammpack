/** @file */

#if defined(__SPAMM_MM_H)
#warn Already included spamm_mm.h
#else

/** Define in case spamm_mm.h has been included. */
#define __SPAMM_MM_H 1

#include "config.h"

/** \page page_mm The SpAMM memory manager.
 *
 * \brief A memory manager to manage contigueous memory allocation.
 *
 * In order to allow for dynamic memory allocation in a contigueous chunk of
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
 *
 * \author Nicolas Bock <nicolasbock@gmail.com>
 * \author Matt Challacombe <matt.challacombe@gmail.com>.
 */

/** A memory chunk.
 */
struct spamm_mm_t
{
  /** Link to previous chunk. */
  struct spamm_mm_t *previous;

  /** Link to next chunk. */
  struct spamm_mm_t *next;

  /** How large are the chunks? */
  unsigned int chunksize;

  /** How many bytes are already allocated in this chunk? */
  unsigned int bytes_allocated;

  /** The data. */
  void *data;
};

struct spamm_mm_t *
spamm_mm_initialize (const unsigned int chunksize);

void *
spamm_mm_allocate (const unsigned int size, struct spamm_mm_t *memory);

#endif

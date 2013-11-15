/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "chunk.h"

#include "index.h"

#include <assert.h>
#include <malloc.h>
#include <stdarg.h>
#include <string.h>

#ifdef DEBUG_OUTPUT
#define DEBUG(message, ...) printf("[%s:%d (%s) DEBUG] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define DEBUG(message, ...) /* stripped DEBUG statement. */
#endif

/** The data layout of a chunk. */
struct chunk_t
{
  /** The size of this matrix block. */
  int blocksize;

  /** The matrix data. */
  double block[0];
};

/** Get the chunksize.
 */
size_t
chunk_sizeof (const int blocksize)
{
  return sizeof(struct chunk_t)+sizeof(double)*blocksize*blocksize;
}

/** Allocate a chunk.
 *
 * @param blocksize The size of the matrix block.
 *
 * @return The newly allocated chunk.
 */
void *
chunk_alloc (const int blocksize)
{
  void *chunk = malloc(chunk_sizeof(blocksize));
  memset(chunk, 0, chunk_sizeof(blocksize));

  struct chunk_t *chunk_ptr = (struct chunk_t*) chunk;
  chunk_ptr->blocksize = blocksize;

  DEBUG("allocating chunk, blocksize = %d, sizeof(chunk) = %lu\n",
      blocksize, chunk_sizeof(blocksize));

  return chunk;
}

/** Set a chunk.
 *
 * @param chunk The chunk.
 * @param A The dense matrix A. This matrix has to have the correct size, i.e.
 * it's blocksize has to match the blocksize of the chunk.
 */
void chunk_set (void *const chunk, const double *const A)
{
  struct chunk_t *chunk_ptr = (struct chunk_t*) chunk;
  memcpy(chunk_ptr->block, A, sizeof(double)*chunk_ptr->blocksize*chunk_ptr->blocksize);
  DEBUG("set chunk, blocksize = %d\n", chunk_ptr->blocksize);
  chunk_print(chunk, "chunk");
}

/** Get the matrix norm of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The square of the Frobenius norm.
 */
double
chunk_get_norm (const void *const chunk)
{
  struct chunk_t *chunk_ptr = (struct chunk_t*) chunk;
  double norm = 0;
  for(int i = 0; i < chunk_ptr->blocksize*chunk_ptr->blocksize; i++)
  {
    norm += chunk_ptr->block[i]*chunk_ptr->block[i];
  }
  DEBUG("chunk norm = %e\n", norm);
  return norm;
}

/** Print a chunk.
 *
 * @param chunk The chunk.
 * @param format The format string. This follows the printf() approach.
 */
void
chunk_print (const void *const chunk,
    const char *const format, ...)
{
  struct chunk_t *chunk_ptr = (struct chunk_t*) chunk;
  char tag[2000];
  va_list ap;

  va_start(ap, format);
  vsnprintf(tag, 2000, format, ap);

  printf("%s\n", tag);
  for(int i = 0; i < chunk_ptr->blocksize; i++) {
    for(int j = 0; j < chunk_ptr->blocksize; j++)
    {
      printf(" % e", chunk_ptr->block[BLOCK_INDEX(i, j, 0, 0, chunk_ptr->blocksize)]);
    }
    printf("\n");
  }
}

/** Add two chunks.
 *
 * @f[ A \leftarrow \alpha A + \beta B @f]
 *
 * @param alpha The scalar alpha.
 * @param A The chunk A.
 * @param beta The scalar beta.
 * @param B The chunk B.
 */
void
chunk_add (const double alpha, void *const A,
    const double beta, const void *const B)
{
  struct chunk_t *c_A = (struct chunk_t*) A;
  struct chunk_t *c_B = (struct chunk_t*) B;

  assert(c_A->blocksize == c_B->blocksize);

  DEBUG("adding two chunks, alpha = %e, beta = %e\n", alpha, beta);

  for(int i = 0; i < c_A->blocksize*c_A->blocksize; i++)
  {
    c_A->block[i] = alpha*c_A->block[i]+beta*c_B->block[i];
  }
}

/** Multiply two chunks.
 *
 * @f[ C \leftarrow A \times B @f]
 *
 * @param A Chunk A.
 * @param B Chunk B.
 * @param C Chunk C.
 */
void
chunk_multiply (const void *const A, const void *const B, void *const C)
{
  struct chunk_t *c_A = (struct chunk_t*) A;
  struct chunk_t *c_B = (struct chunk_t*) B;
  struct chunk_t *c_C = (struct chunk_t*) C;

  assert(c_A->blocksize == c_B->blocksize);
  assert(c_A->blocksize == c_C->blocksize);

  DEBUG("multiplying chunks, blocksize = %d\n", c_A->blocksize);

  for(int i = 0; i < c_A->blocksize; i++) {
    for(int j = 0; j < c_A->blocksize; j++) {
      for(int k = 0; k < c_A->blocksize; k++)
      {
        c_C->block[BLOCK_INDEX(i, j, 0, 0, c_A->blocksize)] +=
          c_A->block[BLOCK_INDEX(i, k, 0, 0, c_A->blocksize)] *
          c_B->block[BLOCK_INDEX(k, j, 0, 0, c_A->blocksize)];
      }
    }
  }
}

/** Get the trace of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The trace of the chunk.
 */
double
chunk_trace (const void *const chunk)
{
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double trace = 0;
  for(int i = 0; i < ptr->blocksize; i++)
  {
    trace += ptr->block[BLOCK_INDEX(i, i, 0, 0, ptr->blocksize)];
  }
  return trace;
}

/** Scale a chunk, i.e. multiply it by a scalar.
 *
 * @f[ A \leftarrow \alpha A @f]
 *
 * @param alpha The scalar alpha.
 * @param A The chunk.
 */
void
chunk_scale (const double alpha, void *const chunk)
{
  DEBUG("FIXME\n");
}

/** Add a scaled identity matrix to a chunk.
 *
 * @f[ A \leftarrow A + \alpha \times \mathrm{Id} @f]
 *
 * @param alpha The scalar alpha.
 * @param chunk The chunk.
 */
void
chunk_add_identity (const double alpha, void *const chunk)
{
  DEBUG("FIXME\n");
}

/** Convert a chunk to a dense matrix.
 *
 * @param chunk The chunk.
 *
 * @return The dense matrix. This matrix is sized to blocksize, which is
 * stored in the chunk.
 */
double *
chunk_to_dense (const void *const chunk)
{
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double *A = calloc(ptr->blocksize*ptr->blocksize, sizeof(double));
  memcpy(A, ptr->block, sizeof(double)*ptr->blocksize*ptr->blocksize);
  return A;
}

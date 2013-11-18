/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "chunk.h"

#include <assert.h>
#include <malloc.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/** A convenience macro for printing some debugging output. */
#ifdef DEBUG_OUTPUT
#define DEBUG(message, ...) printf("[%s:%d (%s) DEBUG] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define DEBUG(message, ...) /* stripped DEBUG statement. */
#endif

/** A convenience macro for printing some info level message. */
#define INFO(message, ...) printf("[%s:%d (%s)] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)

/** A convenience macro for print a fatal error message and terminating the
 * code. */
#define ABORT(message, ...) printf("[%s:%d (%s) FATAL] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__); exit(1)

/** Calculate a row-major offset. */
#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))

/** Calculate a column-major offset. */
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))

/** Calculate the offset into a tiled matrix block. */
#define INDEX(i, j, N) ROW_MAJOR(i, j, N)

/** Convert the norm offset into a pointer. */
#define NORM_POINTER(ptr) (double*) ((intptr_t) \
    ((struct chunk_t*) (ptr))->norm_2 + (intptr_t) ((struct chunk_t*) ptr)->data)

/** Convert the matrix data offset into a pointer. */
#define MATRIX_POINTER(i, j, ptr) (double*) ((intptr_t) \
    ((struct chunk_t*) ptr)->A_block[INDEX(i, j, ((struct chunk_t*) ptr)->N_block)] \
    + (intptr_t) ((struct chunk_t*) ptr)->data)

/** A simple square. */
#define SQUARE(x) (x)*(x)

/** The data layout of a chunk. */
struct chunk_t
{
  /** The lower row bound of this chunk. */
  unsigned int i_lower;

  /** The lower column bound of this chunk. */
  unsigned int j_lower;

  /** The matrix size. The chunk's bounds can extend beyond the matrix size,
   * since we use padding. */
  int N;

  /** The size of this matrix chunk. */
  int N_chunk;

  /** The sub-matrix size. The chunk contains a N_chunk x N_chunk matrix which
   * is tiled in N_basic x N_basic sub-matrices. The submatrices are then
   * multiplied using the SpAMM algorithm, i.e. the norm is used to drop
   * insignificant product contributions at the N_basic level. */
  int N_basic;

  /** The size of the blocked chunk, i.e. N_chunk/N_basic. This is stored here
   * simply for convenience. */
  int N_block;

  /** The N_chunk/N_basic x N_chunk/N_basic size array of norms of the
   * sub-matrices at the N_basic level. The pointer is a relative offset from
   * the start of the data field. */
  double *norm_2;

  /** The N_chunk/N_basic x N_chunk/N_basic size array of pointers to N_basic
   * x N_basic sub-matrices. The pointer is a relative offset from the start
   * of the data field. */
  double **A_block;

  /** The chunk data. This includes the norms and the matrix elements. The
   * exact layout is:
   *
   * N_block*N_block*sizeof(double)
   * N_chunk*N_chunk*sizeof(double)
   * */
  double data[0];
};

/** Get the chunksize.
 *
 * @param N_chunk The size of this matrix chunk.
 * @param N_basic The size of the basic sub-matrices.
 */
size_t
chunk_sizeof (const int N_chunk, const int N_basic)
{
  return sizeof(struct chunk_t)
    +sizeof(double)*N_basic*N_basic
    +sizeof(double)*N_chunk*N_chunk;
}

/** Allocate a chunk.
 *
 * @param N_chunk The size of the matrix chunk.
 * @param N_basic The size of the basic sub-matrices.
 * @param N The matrix size.
 * @param i_lower The lower row bound of this chunk.
 * @param j_lower The lower column bound of this chunk.
 *
 * @return The newly allocated chunk.
 */
void *
chunk_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const unsigned int i_lower,
    const unsigned int j_lower)
{
  void *chunk = malloc(chunk_sizeof(N_chunk, N_basic));
  memset(chunk, 0, chunk_sizeof(N_chunk, N_basic));

  struct chunk_t *ptr = (struct chunk_t*) chunk;

  ptr->i_lower = i_lower;
  ptr->j_lower = j_lower;
  ptr->N = N;
  ptr->N_chunk = N_chunk;
  ptr->N_basic = N_basic;

  if(N_chunk%N_basic != 0)
  {
    ABORT("N_chunk has to be an integer multiply of N_basic\n");
  }

  ptr->N_block = N_chunk/N_basic;

  /* Set the relative offsets. */
  ptr->norm_2 = 0;
  ptr->A_block = (double**) ((intptr_t) ptr->N_block*ptr->N_block*sizeof(double));

  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      ptr->A_block[INDEX(i, j, ptr->N_block)] =
        (double*) (
            (intptr_t) INDEX(i, j, ptr->N_block)
            *ptr->N_basic*ptr->N_basic*sizeof(double));
    }
  }

  DEBUG("allocating chunk at %p, N = %d, N_chunk = %d, "
      "N_basic = %d, sizeof(chunk) = %lu\n",
      chunk, N, N_chunk, N_basic, chunk_sizeof(N_chunk, N_basic));

  return chunk;
}

/** Set a chunk.
 *
 * @param chunk The chunk.
 * @param A The dense matrix A. This matrix has to have the correct size, i.e.
 * this function expects N_chunk x N_chunk elements. The elements of A are
 * expected in column-major order.
 */
void chunk_set (void *const chunk, const double *const A)
{
  struct chunk_t *ptr = (struct chunk_t*) chunk;

  double *norm_2 = NORM_POINTER(chunk);
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      /* Convert offset pointer into actual pointer. */
      double *A_basic = MATRIX_POINTER(i, j, chunk);
      norm_2[INDEX(i, j, ptr->N_block)] = 0;
      for(int k = 0; k < ptr->N_basic; k++) {
        for(int l = 0; l < ptr->N_basic; l++)
        {
          A_basic[INDEX(k, l, ptr->N_basic)] = A[COLUMN_MAJOR(i*ptr->N_block+k, j*ptr->N_block+l, ptr->N_chunk)];
          norm_2[INDEX(i, j, ptr->N_block)] += SQUARE(A_basic[INDEX(k, l, ptr->N_basic)]);
        }
      }
    }
  }

  DEBUG("set chunk at %p, N_chunk = %d\n", chunk, ptr->N_chunk);
#ifdef DEBUG_OUTPUT
  chunk_print(chunk, "chunk");
#endif
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
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double *norm_2 = NORM_POINTER(chunk);
  double norm = 0;
  for(int i = 0; i < SQUARE(ptr->N_block); i++)
  {
    norm += norm_2[i];
  }
  DEBUG("chunk at %p, norm = %e\n", chunk, norm);
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
  struct chunk_t *ptr = (struct chunk_t*) chunk;

  char tag[2000];
  va_list ap;

  va_start(ap, format);
  vsnprintf(tag, 2000, format, ap);

  printf("%s\n", tag);
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = MATRIX_POINTER(i, j, chunk);
      printf("block(%d,%d):\n", i, j);
      for(int k = 0; k < ptr->N_basic; k++) {
        for(int l = 0; l < ptr->N_basic; l++)
        {
          printf(" % e", A_basic[INDEX(k, l, ptr->N_basic)]);
        }
      }
      printf("\n");
    }
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
  struct chunk_t *ptr_A = (struct chunk_t*) A;
  struct chunk_t *ptr_B = (struct chunk_t*) B;

  DEBUG("A at %p, N_chunk = %d, B at %p, N_chunk = %d\n", A, ptr_A->N_chunk, B, ptr_B->N_chunk);

  assert(ptr_A->N_chunk == ptr_B->N_chunk);
  assert(ptr_A->N_basic == ptr_B->N_basic);

  DEBUG("adding two chunks, alpha = %e, beta = %e\n", alpha, beta);

  for(int i = 0; i < ptr_A->N_block; i++) {
    for(int j = 0; j < ptr_A->N_block; j++)
    {
      double *A_basic = MATRIX_POINTER(i, j, A);
      double *B_basic = MATRIX_POINTER(i, j, B);

      for(int k = 0; k < ptr_A->N_basic; k++) {
        for(int l = 0; l < ptr_A->N_basic; l++)
        {
          A_basic[INDEX(k, l, ptr_A->N_basic)] =
            alpha*A_basic[INDEX(k, l, ptr_A->N_basic)]
            +beta*B_basic[INDEX(k, l, ptr_A->N_basic)];
        }
      }
    }
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
  struct chunk_t *ptr_A = (struct chunk_t*) A;
  struct chunk_t *ptr_B = (struct chunk_t*) B;
  struct chunk_t *ptr_C = (struct chunk_t*) C;

  assert(ptr_A->N_chunk == ptr_B->N_chunk);
  assert(ptr_A->N_chunk == ptr_C->N_chunk);

  assert(ptr_A->N_basic == ptr_B->N_basic);
  assert(ptr_A->N_basic == ptr_C->N_basic);

  DEBUG("multiplying chunks A(%p) B(%p) C(%p), N_chunk = %d\n", A, B, C,
      ptr_A->N_chunk);

  /* Simple tiling over the basic sub-matrix blocks. */
  double *norm_2 = NORM_POINTER(C);
  for(int i = 0; i < ptr_A->N_block; i++) {
    for(int j = 0; j < ptr_A->N_block; j++)
    {
      double *C_basic = MATRIX_POINTER(i, j, C);
      memset(C_basic, 0, sizeof(double)*ptr_C->N_basic*ptr_C->N_basic);
      norm_2[INDEX(i, j, ptr_C->N_block)] = 0;
      for(int k = 0; k < ptr_A->N_block; k++)
      {
        /* Multiply the blocks. */
        double *A_basic = MATRIX_POINTER(i, k, A);
        double *B_basic = MATRIX_POINTER(k, j, B);
        for(int i_basic = 0; i_basic < ptr_A->N_basic; i_basic++) {
          for(int j_basic = 0; j_basic < ptr_A->N_basic; j_basic++) {
            for(int k_basic = 0; k_basic < ptr_A->N_basic; k_basic++)
            {
              C_basic[INDEX(i_basic, j_basic, ptr_C->N_basic)] +=
                A_basic[INDEX(i_basic, k_basic, ptr_A->N_basic)]
                *B_basic[INDEX(k_basic, j_basic, ptr_B->N_basic)];
            }
            norm_2[INDEX(i, j, ptr_C->N_block)] += SQUARE(C_basic[INDEX(i_basic, j_basic, ptr_C->N_basic)]);
          }
        }
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
  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = MATRIX_POINTER(i, i, chunk);
    for(int j = 0; j < ptr->N_basic
        && ptr->i_lower+i*ptr->N_basic+j < ptr->N
        && ptr->j_lower+i*ptr->N_basic+j < ptr->N; j++)
    {
      trace += A_basic[INDEX(j, j, ptr->N_basic)];
    }
  }
  DEBUG("calculating trace of chunk at %p, trace = %e\n", chunk, trace);
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
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  DEBUG("scaling block at %p by %e\n", chunk, alpha);
  double *A_chunk = MATRIX_POINTER(0, 0, chunk);
  double *norm_2 = NORM_POINTER(chunk);
  for(int i = 0; i < ptr->N_chunk*ptr->N_chunk; i++)
  {
    A_chunk[i] *= alpha;
  }
  for(int i = 0; i < ptr->N_block; i++)
  {
    norm_2[i] *= SQUARE(alpha);
  }
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
  struct chunk_t *ptr = (struct chunk_t*) chunk;

  assert(ptr->i_lower == ptr->j_lower);

  DEBUG("adding %e Id to chunk at %p, N = %d, i_lower = %d, j_lower = %d\n",
      alpha, chunk, ptr->N, ptr->i_lower, ptr->j_lower);

  double *norm_2 = NORM_POINTER(chunk);
  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = MATRIX_POINTER(i, i, chunk);
    for(int j = 0; j < ptr->N_basic &&
        j+ptr->i_lower < ptr->N &&
        j+ptr->j_lower < ptr->N; j++)
    {
      double old_Aij = A_basic[INDEX(j, j, ptr->N_basic)];
      A_basic[INDEX(j, j, ptr->N_basic)] += alpha;
      norm_2[INDEX(i, i, ptr->N_block)] += SQUARE(A_basic[INDEX(j, j, ptr->N_basic)])
        - SQUARE(old_Aij);
    }
  }
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

  double *A = calloc(ptr->N_chunk*ptr->N_chunk, sizeof(double));
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = MATRIX_POINTER(i, j, ptr);
      for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++) {
        for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
        {
          A[COLUMN_MAJOR(i*ptr->N_block+i_basic, j*ptr->N_block+j_basic, ptr->N_chunk)] =
            A_basic[INDEX(i, j, ptr->N_basic)];
        }
      }
    }
  }
  return A;
}

/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "chunk_block.h"
#include "chunk_tiled.h"

#include <assert.h>
#include <malloc.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** A convenience macro for printing some debugging output. */
#ifdef DEBUG_OUTPUT
#define DEBUG(message, ...) printf("[%s:%d (%s) DEBUG] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define DEBUG(message, ...) /* stripped DEBUG statement. */
#endif

/** A convenience macro for printing some info level message. */
#ifdef _OPENMP
#define INFO(message, ...) printf("[%s:%d (%s) thread %d] " message, __FILE__, __LINE__, __func__, omp_get_thread_num(), ##__VA_ARGS__)
#else
#define INFO(message, ...) printf("[%s:%d (%s)] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#endif

/** A convenience macro for print a fatal error message and terminating the
 * code. */
#define ABORT(message, ...) printf("[%s:%d (%s) FATAL] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__); exit(1)

/** A simple square. */
#define SQUARE(x) ((x)*(x))

/** A cube. */
#define CUBE(x) ((x)*(x)*(x))

/** Calculate a row-major offset. */
#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))

/** Calculate a column-major offset. */
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))

/** The data layout of a chunk. */
struct chunk_tiled_t
{
  /** The total size of this chunk. Used for bounds checks. */
  size_t chunksize;

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

  /** The chunk data. This includes the norms and the matrix elements. The
   * exact layout is:
   *
   * A N_block x N_block array of double norm^2 values: N_block*N_block*sizeof(double)
   * The matrix data: N_chunk*N_chunk*sizeof(double)
   * */
  char data[0];
};

/** Calculate the offset into a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param N The size of the matrix.
 *
 * @return The offset.
 */
static inline size_t
chunk_tiled_matrix_offset (const int i, const int j, const int N)
{
  assert(i >= 0);
  assert(j >= 0);
  assert(i < N);
  assert(j < N);

  size_t index = ROW_MAJOR(i, j, N);
  assert(index < SQUARE(N));
  return index;
}

/** Return a pointer to the start of the norm_2 array.
 *
 * @param chunk The chunk.
 *
 * @return The pointer to the norm_2 array.
 */
static inline double *
chunk_tiled_norm_pointer (const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  return (double*) ((intptr_t) ptr->data);
}

/** Return a pointer to a basic submatrix block.
 *
 * @param i The row pointer.
 * @param j The column pointer.
 * @param chunk The chunk.
 *
 * @return The pointer to the start of a (i,j)th basic submatrix.
 */
static inline double *
chunk_tiled_matrix_pointer (const int i, const int j, const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  assert(i >= 0);
  assert(j >= 0);
  assert(i < ptr->N_block);
  assert(j < ptr->N_block);
  intptr_t offset = (intptr_t) (SQUARE(ptr->N_block)*sizeof(double)) /* The norms. */
    + (intptr_t) (SQUARE(ptr->N_basic)*sizeof(double)
        *chunk_tiled_matrix_offset(i, j, ptr->N_block)); /* Offset to basic submatrix. */
  return (double*) ((intptr_t) ptr->data + offset);
}

/** Get the chunksize.
 *
 * @param N_chunk The size of this matrix chunk.
 * @param N_basic The size of the basic sub-matrices.
 *
 * @return The size of the chunk in bytes.
 */
size_t
chunk_tiled_sizeof (const int N_chunk, const int N_basic)
{
  int N_block = N_chunk/N_basic;

  return sizeof(struct chunk_tiled_t)
    +sizeof(double)*N_block*N_block  /* The norms. */
    +sizeof(double)*N_chunk*N_chunk; /* The matrix elements. */
}

/** Allocate a chunk.
 *
 * The chunk stores a N_chunk x N_chunk dense matrix. The matrix is stored in
 * a matrix of N_basic x N_basic submatrices.
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
chunk_tiled_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const unsigned int i_lower,
    const unsigned int j_lower)
{
  DEBUG("allocating new chunk, N_chunk = %d, N_basic = %d, sizeof(chunk) = %lu\n",
      N_chunk, N_basic, chunk_tiled_sizeof(N_chunk, N_basic));

  void *chunk = calloc(chunk_tiled_sizeof(N_chunk, N_basic), 1);

  DEBUG("setting chunk of size %lu\n", chunk_tiled_sizeof(N_chunk, N_basic));

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;

  ptr->chunksize = chunk_tiled_sizeof(N_chunk, N_basic);
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

  DEBUG("allocating chunk at %p, N = %d, N_chunk = %d, "
      "N_basic = %d, sizeof(chunk) = %lu\n",
      chunk, N, N_chunk, N_basic, chunk_tiled_sizeof(N_chunk, N_basic));

  return chunk;
}

/** Update the norms.
 *
 * @param chunk The chunk.
 */
void
chunk_tiled_set_norm (void *const chunk)
{
  assert(chunk != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  double *A = chunk_tiled_matrix_pointer(0, 0, chunk);
  double *norm = chunk_tiled_norm_pointer(chunk);
#pragma omp parallel for default(none) shared(ptr, norm, A)
  for(int i = 0; i < SQUARE(ptr->N_block); i++)
  {
    norm[i] = 0;
    for(int j = 0; j < SQUARE(ptr->N_basic); j++)
    {
      norm[i] += SQUARE(A[i*SQUARE(ptr->N_basic)+j]);
    }
  }
}

/** Set a chunk.
 *
 * @param chunk The chunk.
 * @param A The dense matrix A. This matrix has to have the correct size, i.e.
 * this function expects N_chunk x N_chunk elements. The elements of A are
 * expected in column-major order.
 */
void
chunk_tiled_set (void *const chunk, const double *const A)
{
  assert(chunk != NULL);
  assert(A != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;

#pragma omp parallel for default(none) shared(ptr)
  for(int index = 0; index < SQUARE(ptr->N_block); index++)
  {
    int j = index;
    int i = j/ptr->N_block;
    j %= ptr->N_block;

    double *A_basic = chunk_tiled_matrix_pointer(i, j, chunk);

    for(int k = 0; k < ptr->N_basic; k++) {
      for(int l = 0; l < ptr->N_basic; l++)
      {
        A_basic[chunk_tiled_matrix_offset(k, l, ptr->N_basic)] =
          A[COLUMN_MAJOR(i*ptr->N_basic+k, j*ptr->N_basic+l, ptr->N_chunk)];
      }
    }
  }

  chunk_tiled_set_norm(chunk);

  DEBUG("set chunk at %p, N_chunk = %d\n", chunk, ptr->N_chunk);
#ifdef PRINT_MATRICES
  chunk_tiled_print(chunk, "chunk");
#endif
}

/** Return the matrix size of the chunk matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_chunk.
 */
int
chunk_tiled_get_N_chunk (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_tiled_t*) chunk)->N_chunk;
}

/** Return the matrix size of the basic matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_basic.
 */
int
chunk_tiled_get_N_basic (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_tiled_t*) chunk)->N_basic;
}

/** Get the square of the matrix norm of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The square of the Frobenius norm.
 */
double
chunk_tiled_get_norm_2 (const void *const chunk)
{
  assert(chunk != NULL);

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  double *norm_2 = chunk_tiled_norm_pointer(chunk);
  double norm = 0;

  for(int i = 0; i < SQUARE(ptr->N_block); i++)
  {
    norm += norm_2[i];
  }
  DEBUG("chunk at %p, norm_2 = %e\n", chunk, norm);
  return norm;
}

/** Print a chunk.
 *
 * @param chunk The chunk.
 * @param format The format string. This follows the printf() approach.
 */
void
chunk_tiled_print (const void *const chunk,
    const char *const format, ...)
{
  assert(chunk != NULL);

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;

  char tag[2000];
  va_list ap;

  va_start(ap, format);
  vsnprintf(tag, 2000, format, ap);

  double *norm = chunk_tiled_norm_pointer(chunk);

  printf("%s\n", tag);
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = chunk_tiled_matrix_pointer(i, j, chunk);
      printf("block(%d,%d), norm = %e:\n", i, j, norm[chunk_tiled_matrix_offset(i, j, ptr->N_block)]);
      double trace = 0;
      for(int k = 0; k < ptr->N_basic; k++)
      {
        trace += A_basic[chunk_tiled_matrix_offset(k, k, ptr->N_basic)];
        for(int l = 0; l < ptr->N_basic; l++)
        {
          printf(" % e", A_basic[chunk_tiled_matrix_offset(k, l, ptr->N_basic)]);
        }
        printf("\n");
      }
      printf("\n");
      printf("block(%d,%d), trace = %e\n", i, j, trace);
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
chunk_tiled_add (const double alpha, void *const A,
    const double beta, const void *const B)
{
  assert(A != NULL);
  assert(B != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr_A = (struct chunk_tiled_t*) A;
  struct chunk_tiled_t *ptr_B = (struct chunk_tiled_t*) B;

  DEBUG("A at %p, N_chunk = %d, B at %p, N_chunk = %d\n", A, ptr_A->N_chunk, B, ptr_B->N_chunk);

  assert(ptr_A->N_chunk == ptr_B->N_chunk);
  assert(ptr_A->N_basic == ptr_B->N_basic);

  DEBUG("adding two chunks, alpha = %e, beta = %e\n", alpha, beta);

  double *A_dense = chunk_tiled_matrix_pointer(0, 0, A);
  double *B_dense = chunk_tiled_matrix_pointer(0, 0, B);

#pragma omp parallel for default(none) shared(ptr_A, A_dense, B_dense)
  for(int i = 0; i < SQUARE(ptr_A->N_chunk); i++)
  {
    A_dense[i] = alpha*A_dense[i]+beta*B_dense[i];
  }

  chunk_tiled_set_norm(A);
}

/** Multiply two chunks using the SpAMM algorithm.
 *
 * @f[ C \leftarrow A \times B @f]
 *
 * @param tolerance The SpAMM tolerance.
 * @param A Chunk A.
 * @param B Chunk B.
 * @param C Chunk C.
 * @param symbolic_only Only go through the symbolic part, i.e. don't multiply the
 * basic blocks. Used only for debugging to get a handle on the performance of
 * the tree.
 */
void
chunk_tiled_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C,
    const short symbolic_only)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr_A = (struct chunk_tiled_t*) A;
  struct chunk_tiled_t *ptr_B = (struct chunk_tiled_t*) B;
  struct chunk_tiled_t *ptr_C = (struct chunk_tiled_t*) C;

  assert(ptr_A->N_chunk == ptr_B->N_chunk);
  assert(ptr_A->N_chunk == ptr_C->N_chunk);

  assert(ptr_A->N_basic == ptr_B->N_basic);
  assert(ptr_A->N_basic == ptr_C->N_basic);

  DEBUG("multiplying chunks A(%p) B(%p) C(%p), N_chunk = %d\n", A, B, C,
      ptr_A->N_chunk);

  INFO("SpAMM tolerance = %e\n", tolerance);

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    {
      INFO("running on %d OpenMP threads\n", omp_get_num_threads());
    }
  }
#else
  INFO("running in serial\n");
#endif

  /* Simple tiling over the basic sub-matrix blocks. */
  int complexity = 0;
  double tolerance_2 = SQUARE(tolerance);
  double *norm_A = chunk_tiled_norm_pointer(A);
  double *norm_B = chunk_tiled_norm_pointer(B);

  /* Reset C. */
  memset(ptr_C->data, 0, sizeof(double)*(SQUARE(ptr_A->N_block)+SQUARE(ptr_A->N_chunk)));

#ifdef _OPENMP
  omp_lock_t *C_lock = malloc(sizeof(omp_lock_t)*SQUARE(ptr_A->N_block));
  for(int i = 0; i < SQUARE(ptr_A->N_block); i++)
  {
    omp_init_lock(&C_lock[i]);
  }
#endif

#pragma omp parallel for default(none) shared(tolerance_2, norm_A, norm_B, ptr_A, ptr_B, ptr_C, C_lock) reduction(+:complexity)
  for(int index = 0; index < CUBE(ptr_A->N_block); index++)
  {
    int i = index;
    int j = i/SQUARE(ptr_A->N_block);
    i %= SQUARE(ptr_A->N_block);
    int k = i/ptr_A->N_block;
    i %= ptr_A->N_block;

    //INFO("index = %d, i = %d, j = %d, k = %d\n", index, i, j, k);

    double *const C_basic = chunk_tiled_matrix_pointer(i, j, C);

    if(norm_A[chunk_tiled_matrix_offset(i, k, ptr_A->N_block)]
        * norm_B[chunk_tiled_matrix_offset(k ,j, ptr_A->N_block)]
        > tolerance_2)
    {
      if(!symbolic_only)
      {
        /* Multiply the blocks. */
#ifdef _OPENMP
        omp_set_lock(&C_lock[chunk_tiled_matrix_offset(i, j, ptr_A->N_block)]);
#endif
        const double *const A_basic = chunk_tiled_matrix_pointer(i, k, A);
        const double *const B_basic = chunk_tiled_matrix_pointer(k, j, B);

        chunk_block_multiply(A_basic, B_basic, C_basic, ptr_A->N_basic);

#ifdef _OPENMP
        omp_unset_lock(&C_lock[chunk_tiled_matrix_offset(i, j, ptr_A->N_block)]);
#endif
      }

      complexity++;
    }
  }

#ifdef _OPENMP
  for(int i = 0; i < SQUARE(ptr_A->N_block); i++)
  {
    omp_destroy_lock(&C_lock[i]);
  }
  free(C_lock);
#endif

  chunk_tiled_set_norm(C);

  INFO("complexity %d out of %d\n", complexity, ptr_A->N_block*ptr_A->N_block*ptr_A->N_block);
}

/** Get the trace of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The trace of the chunk.
 */
double
chunk_tiled_trace (const void *const chunk)
{
  assert(chunk != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  double trace = 0;

#pragma omp parallel for default(none) shared(ptr) reduction(+:trace)
  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = chunk_tiled_matrix_pointer(i, i, chunk);
    for(int j = 0; j < ptr->N_basic
        && ptr->i_lower+i*ptr->N_basic+j < ptr->N
        && ptr->j_lower+i*ptr->N_basic+j < ptr->N; j++)
    {
      trace += A_basic[chunk_tiled_matrix_offset(j, j, ptr->N_basic)];
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
 * @param chunk The chunk.
 */
void
chunk_tiled_scale (const double alpha, void *const chunk)
{
  assert(chunk != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  DEBUG("scaling block at %p by %e\n", chunk, alpha);

  double *A_dense = chunk_tiled_matrix_pointer(0, 0, chunk);
#pragma omp parallel for default(none) shared(ptr, A_dense)
  for(int i = 0; i < SQUARE(ptr->N_chunk); i++)
  {
    A_dense[i] *= alpha;
  }

  double *norm_2 = chunk_tiled_norm_pointer(chunk);
#pragma omp parallel for default(none) shared(ptr, norm_2)
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
chunk_tiled_add_identity (const double alpha, void *const chunk)
{
  assert(chunk != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;

  assert(ptr->i_lower == ptr->j_lower);

  DEBUG("adding %e Id to chunk at %p, N = %d, i_lower = %d, j_lower = %d\n",
      alpha, chunk, ptr->N, ptr->i_lower, ptr->j_lower);

#pragma omp parallel for default(none) shared(ptr)
  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = chunk_tiled_matrix_pointer(i, i, chunk);
    for(int j = 0; j < ptr->N_basic &&
        j+ptr->i_lower < ptr->N &&
        j+ptr->j_lower < ptr->N; j++)
    {
      A_basic[chunk_tiled_matrix_offset(j, j, ptr->N_basic)] += alpha;
    }
  }

  chunk_tiled_set_norm(chunk);
}

/** Convert a chunk to a dense matrix.
 *
 * @param chunk The chunk.
 *
 * @return The dense matrix. This matrix is sized to blocksize, which is
 * stored in the chunk.
 */
double *
chunk_tiled_to_dense (const void *const chunk)
{
  assert(chunk != NULL);

#ifdef FORCE_OMP_NUM_THREADS
  assert(omp_get_max_threads() == FORCE_OMP_NUM_THREADS);
#endif

  struct chunk_tiled_t *ptr = (struct chunk_tiled_t*) chunk;
  double *A = malloc(SQUARE(ptr->N_chunk)*sizeof(double));

#pragma omp parallel for default(none) shared(ptr, A)
  for(int index = 0; index < SQUARE(ptr->N_block); index++)
  {
    int j = index;
    int i = j/ptr->N_block;
    j %= ptr->N_block;

    double *A_basic = chunk_tiled_matrix_pointer(i, j, ptr);
    for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++) {
      for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
      {
        A[COLUMN_MAJOR(i*ptr->N_basic+i_basic, j*ptr->N_basic+j_basic, ptr->N_chunk)] =
          A_basic[chunk_tiled_matrix_offset(i_basic, j_basic, ptr->N_basic)];
      }
    }
  }

  return A;
}

/** Delete the chunk.
 *
 * @param chunk The chunk.
 */
void
chunk_tiled_delete (void **const chunk)
{
  free(*chunk);
  *chunk = NULL;
}

/** Get the complexity of the last chunk operation.
 *
 * The complexity is simply the number of block operations performed. For the
 * multiply this is in the worst case \f$ \left( N_{chunk}/N_{basic}
 * \right)^{3} \f$, but is smaller for matrices with decay and a tolerance \f$
 * \tau > 0 \f$.
 *
 * @param chunk The chunk.
 *
 * @return The complexity count.
 */
size_t
chunk_tiled_get_complexity (const void *const chunk)
{
  ABORT("FIXME\n");
}

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
#define INFO(message, ...) printf("[%s:%d (%s)] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)

/** A convenience macro for print a fatal error message and terminating the
 * code. */
#define ABORT(message, ...) printf("[%s:%d (%s) FATAL] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__); exit(1)

/** A simple square. */
#define SQUARE(x) (x)*(x)

/** A cube. */
#define CUBE(x) (x)*(x)*(x)

/** Calculate a row-major offset. */
#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))

/** Calculate a column-major offset. */
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))

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

  /** The chunk data. This includes the norms and the matrix elements. The
   * exact layout is:
   *
   * A N_block x N_block array of double norm^2 values: N_block*N_block*sizeof(double)
   * The matrix data: N_chunk*N_chunk*sizeof(double)
   * */
  char data[0];
};

struct chunk_linear_state_t
{
  size_t index;
  int i;
  int j;
  int k;
};

void
chunk_init_linear_state (struct chunk_linear_state_t *const state)
{
  INFO("initializing state at %p\n", state);
  memset(state, 0, sizeof(struct chunk_linear_state_t));
}

/** Map a linear index in 3-D convolution space to three indices.
 *
 * @param index The linear index.
 * @param i The index i.
 * @param j The index j.
 * @param k The index k.
 * @param N The maximum index in each direction.
 * @param state The state.
 */
void
chunk_map_linear_index (const size_t index,
    int *const i,
    int *const j,
    int *const k,
    const int N,
    struct chunk_linear_state_t *state)
{
  int bitmask = 1;

  *i = state->i;
  *j = state->j;
  *k = state->k;

  for( ; state->index <= index; state->index++)
  {
    if(state->index > 0)
    {
      *i += 1;
      if((*i)/N > 0)
      {
        *j += (*i)/N;
        *i = (*i)%N;

        if((*j)/N > 0)
        {
          *k += (*j)/N;
          *j = (*j)%N;

          if((*k)/N > 0)
          {
            ABORT("error\n");
          }
        }
      }
    }
  }

  state->i = *i;
  state->j = *j;
  state->k = *k;

  INFO("state at %p, linear index = %lu, index = { %d, %d, %d }, N = %d\n",
      state, index, *i, *j, *k, N);
}

/** Calculate the offset into a tiled matrix block. */
size_t
chunk_index (const int i, const int j, const int N)
{
  assert(i < N);
  assert(j < N);

  size_t index = ROW_MAJOR(i, j, N);
  assert(index < SQUARE(N));
  return index;
}

/** Return a pointer to the start of the norm_2 array. */
double *
chunk_norm_pointer (const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  return (double*) ((intptr_t) ptr->data);
}

/** Return a pointer to a basic submatrix block. */
double *
chunk_matrix_pointer (const int i, const int j, const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_t *ptr = (struct chunk_t*) chunk;
  return (double*) ((intptr_t) ptr->data
      + (intptr_t) SQUARE(ptr->N_block)*sizeof(double)
      + (intptr_t) SQUARE(ptr->N_basic)*sizeof(double)*chunk_index(i, j, ptr->N_block));
}

/** Get the chunksize.
 *
 * @param N_chunk The size of this matrix chunk.
 * @param N_basic The size of the basic sub-matrices.
 */
size_t
chunk_sizeof (const int N_chunk, const int N_basic)
{
  int N_block = N_chunk/N_basic;

  return sizeof(struct chunk_t)
    +sizeof(double)*N_block*N_block  /* The norms. */
    +sizeof(double)*N_chunk*N_chunk; /* The matrix elements. */
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
  DEBUG("allocating new chunk, N_chunk = %d, N_basic = %d, sizeof(chunk) = %lu\n",
      N_chunk, N_basic, chunk_sizeof(N_chunk, N_basic));

  void *chunk = calloc(chunk_sizeof(N_chunk, N_basic), 1);

  DEBUG("setting chunk of size %lu\n", chunk_sizeof(N_chunk, N_basic));

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
  assert(chunk != NULL);
  assert(A != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;

  double *norm_2 = chunk_norm_pointer(chunk);
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = chunk_matrix_pointer(i, j, chunk);

      int index_norm = chunk_index(i, j, ptr->N_block);
      assert(index_norm < SQUARE(ptr->N_block));

      norm_2[index_norm] = 0;
      for(int k = 0; k < ptr->N_basic; k++) {
        for(int l = 0; l < ptr->N_basic; l++)
        {
          int index_basic = chunk_index(k, l, ptr->N_basic);
          int index_dense = COLUMN_MAJOR(i*ptr->N_basic+k, j*ptr->N_basic+l, ptr->N_chunk);

          assert(index_basic < SQUARE(ptr->N_basic));
          assert(index_dense < SQUARE(ptr->N_chunk));

          A_basic[index_basic] = A[index_dense];
          norm_2[index_norm] += SQUARE(A_basic[index_basic]);
        }
      }
    }
  }

  DEBUG("set chunk at %p, N_chunk = %d\n", chunk, ptr->N_chunk);
#ifdef PRINT_MATRICES
  chunk_print(chunk, "chunk");
#endif
}

/** Return the matrix size of the chunk matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_chunk.
 */
int
chunk_get_N_chunk (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_t*) chunk)->N_chunk;
}

/** Return the matrix size of the basic matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_basic.
 */
int
chunk_get_N_basic (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_t*) chunk)->N_basic;
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double *norm_2 = chunk_norm_pointer(chunk);
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;

  char tag[2000];
  va_list ap;

  va_start(ap, format);
  vsnprintf(tag, 2000, format, ap);

  printf("%s\n", tag);
  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = chunk_matrix_pointer(i, j, chunk);
      printf("block(%d,%d):\n", i, j);
      for(int k = 0; k < ptr->N_basic; k++) {
        for(int l = 0; l < ptr->N_basic; l++)
        {
          printf(" % e", A_basic[chunk_index(k, l, ptr->N_basic)]);
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
  assert(A != NULL);
  assert(B != NULL);

  struct chunk_t *ptr_A = (struct chunk_t*) A;
  struct chunk_t *ptr_B = (struct chunk_t*) B;

  DEBUG("A at %p, N_chunk = %d, B at %p, N_chunk = %d\n", A, ptr_A->N_chunk, B, ptr_B->N_chunk);

  assert(ptr_A->N_chunk == ptr_B->N_chunk);
  assert(ptr_A->N_basic == ptr_B->N_basic);

  DEBUG("adding two chunks, alpha = %e, beta = %e\n", alpha, beta);

  for(int i = 0; i < ptr_A->N_block; i++) {
    for(int j = 0; j < ptr_A->N_block; j++)
    {
      double *A_basic = chunk_matrix_pointer(i, j, A);
      double *B_basic = chunk_matrix_pointer(i, j, B);

      for(int k = 0; k < ptr_A->N_basic; k++) {
        for(int l = 0; l < ptr_A->N_basic; l++)
        {
          A_basic[chunk_index(k, l, ptr_A->N_basic)] =
            alpha*A_basic[chunk_index(k, l, ptr_A->N_basic)]
            +beta*B_basic[chunk_index(k, l, ptr_A->N_basic)];
        }
      }
    }
  }
}

/** Multiply two chunks using the SpAMM algorithm.
 *
 * @f[ C \leftarrow A \times B @f]
 *
 * @param tolerance The SpAMM tolerance.
 * @param A Chunk A.
 * @param B Chunk B.
 * @param C Chunk C.
 */
void
chunk_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

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
  int complexity = 0;
  double tolerance_2 = SQUARE(tolerance);
  double *norm_A = chunk_norm_pointer(A);
  double *norm_B = chunk_norm_pointer(B);
  double *norm_C = chunk_norm_pointer(C);

#ifdef _OPENMP
  omp_lock_t C_lock[SQUARE(ptr_A->N_block)];
  for(int i = 0; i < SQUARE(ptr_A->N_block); i++)
  {
    omp_init_lock(&C_lock[i]);
  }
#endif

  short initialized_state = 0;

#pragma omp parallel for default(none) shared(tolerance_2, norm_A, norm_B, norm_C, ptr_A, ptr_B, ptr_C, C_lock) private(initialized_state) reduction(+:complexity)
  for(size_t index = 0; index < CUBE(ptr_A->N_block); index++)
  {
    int i = 0;
    int j = 0;
    int k = 0;

    struct chunk_linear_state_t state;

    if(initialized_state == 0)
    {
      chunk_init_linear_state(&state);
      initialized_state = 1;
    }

    chunk_map_linear_index(index, &i, &j, &k, ptr_A->N_block, &state);

    double *C_basic = chunk_matrix_pointer(i, j, C);
    if(k == 0)
    {
      memset(C_basic, 0, sizeof(double)*ptr_C->N_basic*ptr_C->N_basic);
      norm_C[chunk_index(i, j, ptr_C->N_block)] = 0;
    }

    if(norm_A[chunk_index(i, k, ptr_A->N_block)]
        * norm_B[chunk_index(k ,j, ptr_A->N_block)]
        > tolerance_2)
    {
      /* Multiply the blocks. */
#ifdef _OPENMP
      omp_set_lock(&C_lock[chunk_index(i, j, ptr_A->N_block)]);
#endif
      double *A_basic = chunk_matrix_pointer(i, k, A);
      double *B_basic = chunk_matrix_pointer(k, j, B);
      for(int i_basic = 0; i_basic < ptr_A->N_basic; i_basic++) {
        for(int j_basic = 0; j_basic < ptr_A->N_basic; j_basic++) {
          for(int k_basic = 0; k_basic < ptr_A->N_basic; k_basic++)
          {
            C_basic[chunk_index(i_basic, j_basic, ptr_C->N_basic)] +=
              A_basic[chunk_index(i_basic, k_basic, ptr_A->N_basic)]
              *B_basic[chunk_index(k_basic, j_basic, ptr_B->N_basic)];
          }
          norm_C[chunk_index(i, j, ptr_C->N_block)] += SQUARE(C_basic[chunk_index(i_basic, j_basic, ptr_C->N_basic)]);
        }
      }
#ifdef _OPENMP
      omp_unset_lock(&C_lock[chunk_index(i, j, ptr_A->N_block)]);
#endif

      complexity++;
    }
  }

  DEBUG("complexity %d out of %d\n", complexity, ptr_A->N_block*ptr_A->N_block*ptr_A->N_block);
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double trace = 0;

  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = chunk_matrix_pointer(i, i, chunk);
    for(int j = 0; j < ptr->N_basic
        && ptr->i_lower+i*ptr->N_basic+j < ptr->N
        && ptr->j_lower+i*ptr->N_basic+j < ptr->N; j++)
    {
      trace += A_basic[chunk_index(j, j, ptr->N_basic)];
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;
  DEBUG("scaling block at %p by %e\n", chunk, alpha);
  double *A_chunk = chunk_matrix_pointer(0, 0, chunk);
  double *norm_2 = chunk_norm_pointer(chunk);
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;

  assert(ptr->i_lower == ptr->j_lower);

  DEBUG("adding %e Id to chunk at %p, N = %d, i_lower = %d, j_lower = %d\n",
      alpha, chunk, ptr->N, ptr->i_lower, ptr->j_lower);

  double *norm_2 = chunk_norm_pointer(chunk);
  for(int i = 0; i < ptr->N_block; i++)
  {
    double *A_basic = chunk_matrix_pointer(i, i, chunk);
    for(int j = 0; j < ptr->N_basic &&
        j+ptr->i_lower < ptr->N &&
        j+ptr->j_lower < ptr->N; j++)
    {
      double old_Aij = A_basic[chunk_index(j, j, ptr->N_basic)];
      A_basic[chunk_index(j, j, ptr->N_basic)] += alpha;
      norm_2[chunk_index(i, i, ptr->N_block)] += SQUARE(A_basic[chunk_index(j, j, ptr->N_basic)])
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
  assert(chunk != NULL);

  struct chunk_t *ptr = (struct chunk_t*) chunk;
  double *A = calloc(ptr->N_chunk*ptr->N_chunk, sizeof(double));

  for(int i = 0; i < ptr->N_block; i++) {
    for(int j = 0; j < ptr->N_block; j++)
    {
      double *A_basic = chunk_matrix_pointer(i, j, ptr);
      for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++) {
        for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
        {
          A[COLUMN_MAJOR(i*ptr->N_basic+i_basic, j*ptr->N_basic+j_basic, ptr->N_chunk)] =
            A_basic[chunk_index(i_basic, j_basic, ptr->N_basic)];
        }
      }
    }
  }

  return A;
}

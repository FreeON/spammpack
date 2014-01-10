/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** A convenience macro for printing some debugging output. */
#ifdef DEBUG_OUTPUT
#ifdef _OPENMP
#define DEBUG(message, ...) printf("[%s:%d (%s) thread %d DEBUG] " message, __FILE__, __LINE__, __func__, omp_get_thread_num(), ##__VA_ARGS__)
#else
#define DEBUG(message, ...) printf("[%s:%d (%s) DEBUG] " message, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#endif
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

struct chunk_tree_t
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

  /** The tree depth. The tree starts at tier 0, and grows to 4^{depth}
   * submatrices. */
  int depth;

  /** The offset from the start of the chunk to the array of sub-matrices. */
  intptr_t submatrix_offset;

  /** The chunk data. This includes the norms and the matrix elements. The
   * exact layout is:
   *
   * The tree nodes:
   * tier = 0: root
   * tier = 1: node[4]
   * ...
   * tier = depth: node[4^{depth}]
   * An array of N_basic x N_basic submatrices, the actual matrix data.
   * */
  char data[0];
};

struct chunk_tree_node_t
{
  /** The square of the Frobenius norm of this node. */
  double norm_2;

  /** The size of the basic submatrix. */
  int N_basic;

#ifdef _OPENMP
  /** The OpenMP lock on the matrix. */
  omp_lock_t matrix_lock;
#endif

  union
  {
    /** The children nodes. */
    struct chunk_tree_node_t *child[4];

    /** The submatrix data for leaf nodes. */
    double *matrix;
  }
  data;
};

/** A convenience function to get 2^i, where i is integer.
 *
 * @param i The integer.
 *
 * @return 2 raised to the power of i.
 */
int
ipow2 (const int i)
{
  if(i == 0)
  {
    return 1;
  }

  else if(i > 0)
  {
    return (1 << i);
  }

  else
  {
    ABORT("argument has to be >= 0\n");
  }
}

/** Get the tree depth.
 *
 * @param N_chunk The size of the matrix stored in the chunk.
 * @param N_basic The basic sub-matrices.
 *
 * @return The depth.
 */
int
chunk_tree_get_depth (const int N_chunk, const int N_basic)
{
  int depth = (int) ceil(log(N_chunk/(double) N_basic)/log(2.));
  if(N_basic*ipow2(depth) != N_chunk)
  {
    ABORT("logic error: N_chunk = %d, N_basic = %d, "
        "depth = %d, N_basic*2^depth = %d\n",
        N_chunk, N_basic, depth, N_basic*ipow2(depth));
  }
  return depth;
}

/** Get the chunksize.
 *
 * @param N_chunk The size of this matrix chunk.
 * @param N_basic The size of the basic sub-matrices.
 *
 * @return The size of the chunk in bytes.
 */
size_t
chunk_tree_sizeof (const int N_chunk, const int N_basic)
{
  size_t matrix_size = SQUARE(N_chunk)*sizeof(double);
  size_t tree_size = 0;
  for(int tier = 0; tier <= chunk_tree_get_depth(N_chunk, N_basic); tier++)
  {
    tree_size += ipow2(2*tier)*sizeof(struct chunk_tree_node_t);
  }

  return sizeof(struct chunk_tree_t) + matrix_size + tree_size;
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
chunk_tree_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const int i_lower,
    const int j_lower)
{
  void *chunk = calloc(chunk_tree_sizeof(N_chunk, N_basic), 1);
  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;

  ptr->chunksize = chunk_tree_sizeof(N_chunk, N_basic);
  ptr->i_lower = i_lower;
  ptr->j_lower = j_lower;
  ptr->N = N;
  ptr->N_chunk = N_chunk;
  ptr->N_basic = N_basic;
  ptr->depth = chunk_tree_get_depth(N_chunk, N_basic);

  DEBUG("ptr = %p\n", ptr);
  DEBUG("ptr->data = %p\n", ptr->data);
  DEBUG("ptr->data + chunksize = %p\n", (void*) ((intptr_t) ptr + ptr->chunksize));
  DEBUG("sizeof(struct chunk_tree_node_t) = 0x%lx\n", sizeof(struct chunk_tree_node_t));
  DEBUG("sizeof(submatrix) = 0x%lx\n", SQUARE(N_basic)*sizeof(double));

  /* Fill the chunk with a matrix tree. */
  DEBUG("linking tree nodes...\n");
  intptr_t tier_ptr = (intptr_t) ptr->data;
  for(int tier = 0; tier < ptr->depth; tier++)
  {
    intptr_t next_tier_ptr = tier_ptr + ipow2(2*tier)*sizeof(struct chunk_tree_node_t);

    for(int i = 0; i < ipow2(tier); i++)
    {
      for(int j = 0; j < ipow2(tier); j++)
      {
        size_t offset = ROW_MAJOR(i, j, ipow2(tier))*sizeof(struct chunk_tree_node_t);
        struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) (tier_ptr + offset);

        /* Initialize node. */
        node->N_basic = N_basic;

        DEBUG("%d:node(%d,%d) at %p\n", tier, i, j, node);

        for(int i_child = 0; i_child < 2; i_child++)
        {
          for(int j_child = 0; j_child < 2; j_child++)
          {
            size_t child_offset = ROW_MAJOR(2*i+i_child, 2*j+j_child, ipow2(tier+1))*sizeof(struct chunk_tree_node_t);
            struct chunk_tree_node_t *child = (struct chunk_tree_node_t*) (next_tier_ptr + child_offset);
            node->data.child[ROW_MAJOR(i_child, j_child, 2)] = child;
            DEBUG("%d:child(%d,%d) -> %d:node(%d,%d) at %p\n",
                tier, i_child, j_child, tier+1, 2*i+i_child, 2*j+j_child,
                child);
          }
        }
      }
    }

    tier_ptr = next_tier_ptr;
  }

  /* Store pointers to the basic sub-matrices. */
  DEBUG("linking submatrices...\n");
  intptr_t submatrix_ptr = tier_ptr + ipow2(2*ptr->depth)*sizeof(struct chunk_tree_node_t);
  ptr->submatrix_offset = submatrix_ptr - (intptr_t) ptr;
  DEBUG("submatrix_offset = %ld\n", ptr->submatrix_offset);
  for(int i = 0; i < ipow2(ptr->depth); i++)
  {
    for(int j = 0; j < ipow2(ptr->depth); j++)
    {
      size_t offset = ROW_MAJOR(i, j, ipow2(ptr->depth))*sizeof(struct chunk_tree_node_t);
      size_t matrix_offset = ROW_MAJOR(i, j, ipow2(ptr->depth))*SQUARE(N_basic)*sizeof(double);
      struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) (tier_ptr + offset);

      /* Initialize node. */
      node->N_basic = N_basic;
#ifdef _OPENMP
      omp_init_lock(&node->matrix_lock);
#endif

      node->data.matrix = (double*) (submatrix_ptr + matrix_offset);
      DEBUG("%d:node(%d,%d) at %p, submatrix at %p\n", ptr->depth, i, j, node, &(*node->data.matrix));
    }
  }

  DEBUG("chunk ends at %p, done\n", (void *) (submatrix_ptr
        + ipow2(2*ptr->depth)*SQUARE(N_basic)*sizeof(double)));

  return chunk;
}

/** Get the pointer to a sub-matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param chunk The chunk.
 *
 * @return The pointer to the sub-matrix.
 */
double *
chunk_tree_get_submatrix (const int i,
    const int j,
    const struct chunk_tree_t *const chunk)
{
  assert(chunk != NULL);

  assert(i >= 0);
  assert(j >= 0);

  assert(i < chunk->N_chunk/chunk->N_basic);
  assert(j < chunk->N_chunk/chunk->N_basic);

  double *submatrix_ptr = (double*) ((intptr_t) chunk
      + chunk->submatrix_offset
      + ROW_MAJOR(i, j, ipow2(chunk->depth))*SQUARE(chunk->N_basic)*sizeof(double));

  return submatrix_ptr;
}

/** Get the pointer to a tree node in the chunk.
 *
 * @param tier The tier.
 * @param i The row index.
 * @param j The column index.
 * @param chunk The chunk.
 *
 * @return The pointer to the tree node.
 */
struct chunk_tree_node_t *
chunk_tree_get_node (const int tier,
    const int i,
    const int j,
    const struct chunk_tree_t *const chunk)
{
  assert(chunk != NULL);
  assert(tier >= 0);
  assert(tier <= chunk->depth);
  assert(i >= 0);
  assert(j >= 0);
  assert(i < ipow2(tier));
  assert(j < ipow2(tier));

  intptr_t offset = 0;
  for(int i_tier = 0; i_tier < tier; i_tier++)
  {
    offset += ipow2(2*i_tier)*sizeof(struct chunk_tree_node_t);
  }
  offset += ROW_MAJOR(i, j, ipow2(tier))*sizeof(struct chunk_tree_node_t);
  struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) ((intptr_t) chunk->data + offset);
  DEBUG("%d:node(%d,%d) at %p\n", tier, i, j, node);
  return node;
}

/** Update the norms in the tree part of the chunk.
 *
 * @param chunk The chunk.
 */
void
chunk_tree_update_norm (struct chunk_tree_t *const chunk)
{
  assert(chunk != NULL);

  DEBUG("depth = %d\n", chunk->depth);
  for(int tier = chunk->depth; tier >= 0; tier--)
  {
    for(int i = 0; i < ipow2(tier); i++)
    {
      for(int j = 0; j < ipow2(tier); j++)
      {
        struct chunk_tree_node_t *node = chunk_tree_get_node(tier, i, j, chunk);

        if(tier == chunk->depth)
        {
          double *submatrix = chunk_tree_get_submatrix(i, j, chunk);

          for(int i_basic = 0; i_basic < chunk->N_basic; i_basic++)
          {
            for(int j_basic = 0; j_basic < chunk->N_basic; j_basic++)
            {
              node->norm_2 += SQUARE(submatrix[ROW_MAJOR(i_basic, j_basic, chunk->N_basic)]);
            }
          }
        }

        else
        {
          node->norm_2 = 0;
          for(int i_child = 0; i_child < 2; i_child++)
          {
            for(int j_child = 0; j_child < 2; j_child++)
            {
              struct chunk_tree_node_t *child = chunk_tree_get_node(tier+1, 2*i+i_child, 2*j+j_child, chunk);
              node->norm_2 += child->norm_2;
            }
          }
        }

        DEBUG("%d:node(%d,%d), norm_2 = %e\n", tier, i, j, node->norm_2);
      }
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
chunk_tree_set (void *const chunk, const double *const A)
{
  struct chunk_tree_t *ptr =  (struct chunk_tree_t*) chunk;

  DEBUG("setting chunk\n");

  /* Store the matrix elements. The input matrix is colum-major ordered. The
   * matrix in the chunk is stored in basic submatrix blocks. We need to copy
   * the matrix elements.
   */
  for(int i = 0; i < ptr->N_chunk/ptr->N_basic; i++)
  {
    for(int j = 0; j < ptr->N_chunk/ptr->N_basic; j++)
    {
      double *submatrix_ptr = chunk_tree_get_submatrix(i, j, chunk);

      for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++)
      {
        for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
        {
          submatrix_ptr[ROW_MAJOR(i_basic, j_basic, ptr->N_basic)] =
            A[COLUMN_MAJOR(i*ptr->N_basic+i_basic, j*ptr->N_basic+j_basic, ptr->N_chunk)];
        }
      }
    }
  }

#ifdef DEBUG_OUTPUT
  DEBUG("A:\n");
  for(int i = 0; i < ptr->N_chunk; i++)
  {
    DEBUG("");
    for(int j = 0; j < ptr->N_chunk; j++)
    {
      printf(" % 1.3f", A[COLUMN_MAJOR(i, j, ptr->N_chunk)]);
    }
    printf("\n");
  }

  DEBUG("sub-matrices:\n");
  for(int i = 0; i < ptr->N_chunk/ptr->N_basic; i++)
  {
    for(int j = 0; j < ptr->N_chunk/ptr->N_basic; j++)
    {
      double *submatrix_ptr = chunk_tree_get_submatrix(i, j, chunk);

      DEBUG("A(%d,%d) at %p\n", i, j, submatrix_ptr);
      for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++)
      {
        DEBUG("");
        for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
        {
          printf(" % 1.3f", submatrix_ptr[ROW_MAJOR(i_basic, j_basic, ptr->N_basic)]);
        }
        printf("\n");
      }
    }
  }
#endif

  /** Update norms. */
  chunk_tree_update_norm(chunk);
}

/** Set a chunk to zero.
 *
 * @param chunk The chunk.
 */
void
chunk_tree_set_zero(struct chunk_tree_t *const chunk)
{
  assert(chunk != NULL);

  DEBUG("setting chunk to zero\n");
  double *matrix = chunk_tree_get_submatrix(0, 0, chunk);
  memset(matrix, 0, SQUARE(chunk->N_chunk)*sizeof(double));
  chunk_tree_update_norm(chunk);
}

/** Get the matrix norm of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The square of the Frobenius norm.
 */
double
chunk_tree_get_norm (const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;
  struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;
  return root->norm_2;
}

/** Multiply two tree nodes using the SpAMM algorithm.
 *
 * @f[ C \leftarrow A \times B @f]
 *
 * @param tolerance The SpAMM tolerance.
 * @param tier The tier.
 * @param depth The tree depth.
 * @param A Chunk A.
 * @param B Chunk B.
 * @param C Chunk C.
 * @param tree_only Only go through the symbolic part, i.e. don't multiply the
 * basic blocks. Used only for debugging to get a handle on the performance of
 * the tree.
 */
void
chunk_tree_multiply_node (const double tolerance_2,
    const int tier,
    const int depth,
    const struct chunk_tree_node_t *const A,
    const struct chunk_tree_node_t *const B,
    struct chunk_tree_node_t *const C,
    const short tree_only)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  DEBUG("%d(%d) A at %p, B at %p, C at %p\n", tier, depth, A, B, C);

  DEBUG("tree_only = %d\n", tree_only);

  if(tier == depth)
  {
    DEBUG("%d(%d) A at %p, B at %p, C at %p\n", tier, depth, A, B, C);

    if(!tree_only)
    {
      double *A_submatrix = A->data.matrix;
      double *B_submatrix = B->data.matrix;
      double *C_submatrix = C->data.matrix;

      DEBUG("mulitplying submatrix, A at %p, B at %p, C at %p\n", &(*A_submatrix), &(*B_submatrix), &(*C_submatrix));

#ifdef _OPENMP
      omp_set_lock(&C->matrix_lock);
#endif

      for(int i = 0; i < A->N_basic; i++)
      {
        for(int j = 0; j < A->N_basic; j++)
        {
          for(int k = 0; k < A->N_basic; k++)
          {
            C_submatrix[ROW_MAJOR(i, j, A->N_basic)] +=
              A_submatrix[ROW_MAJOR(i, k, A->N_basic)]
              * B_submatrix[ROW_MAJOR(k, j, A->N_basic)];
          }
        }
      }

#ifdef _OPENMP
      omp_unset_lock(&C->matrix_lock);
#endif
      DEBUG("done multiplying\n");
    }
  }

  else
  {
    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        for(int k = 0; k < 2; k++)
        {
#pragma omp task default(none) firstprivate(i, j, k) untied
          {
            DEBUG("(%d,%d,%d) starting new task\n", i, j, k);

            struct chunk_tree_node_t *A_child = A->data.child[ROW_MAJOR(i, k, 2)];
            struct chunk_tree_node_t *B_child = B->data.child[ROW_MAJOR(k, j, 2)];
            struct chunk_tree_node_t *C_child = C->data.child[ROW_MAJOR(i, j, 2)];

            DEBUG("(%d,%d,%d) A at %p, B at %p, C at %p\n", i, j, k, A_child, B_child, C_child);

            if(A_child->norm_2*B_child->norm_2 > tolerance_2)
            {
              chunk_tree_multiply_node(tolerance_2, tier+1, depth, A_child, B_child, C_child, tree_only);
            }
          }
        }
      }
    }
    DEBUG("done submitting all tasks, waiting...\n");
#pragma omp taskwait
    DEBUG("done\n");
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
 * @param tree_only Only go through the symbolic part, i.e. don't multiply the
 * basic blocks. Used only for debugging to get a handle on the performance of
 * the tree.
 */
void
chunk_tree_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C,
    const short tree_only)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  struct chunk_tree_t *A_ptr = (struct chunk_tree_t*) A;
  struct chunk_tree_t *B_ptr = (struct chunk_tree_t*) B;
  struct chunk_tree_t *C_ptr = (struct chunk_tree_t*) C;

  /* Reset C. */
  chunk_tree_set_zero(C_ptr);

  double tolerance_2 = SQUARE(tolerance);

  struct chunk_tree_node_t *A_root = (struct chunk_tree_node_t*) A_ptr->data;
  struct chunk_tree_node_t *B_root = (struct chunk_tree_node_t*) B_ptr->data;
  struct chunk_tree_node_t *C_root = (struct chunk_tree_node_t*) C_ptr->data;

  INFO("tree_only = %d\n", tree_only);

  INFO("%dx%d blocked matrix, potentially %d products to consider\n",
      ipow2(A_ptr->depth), ipow2(A_ptr->depth), CUBE(ipow2(A_ptr->depth)));

  if(A_root->norm_2*B_root->norm_2 > tolerance_2)
  {
#pragma omp parallel
    {
#pragma omp master
      {
#ifdef _OPENMP
        INFO("running on %d OpenMP threads\n", omp_get_num_threads());
#endif
#pragma omp task untied
        {
          chunk_tree_multiply_node(tolerance_2, 0, A_ptr->depth, A_root, B_root, C_root, tree_only);
        }
#pragma omp taskwait
      }
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
chunk_tree_to_dense (const void *const chunk)
{
  assert(chunk != NULL);

  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;

  double *A = calloc(SQUARE(ptr->N_chunk), sizeof(double));
  for(int i = 0; i < ptr->N_chunk/ptr->N_basic; i++)
  {
    for(int j = 0; j < ptr->N_chunk/ptr->N_basic; j++)
    {
      double *submatrix = chunk_tree_get_submatrix(i, j, chunk);

      for(int i_basic = 0; i_basic < ptr->N_basic; i_basic++)
      {
        for(int j_basic = 0; j_basic < ptr->N_basic; j_basic++)
        {
          A[COLUMN_MAJOR(i*ptr->N_basic+i_basic, j*ptr->N_basic+j_basic, ptr->N_chunk)] =
            submatrix[ROW_MAJOR(i_basic, j_basic, ptr->N_basic)];
        }
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
chunk_tree_delete (void **const chunk)
{
  ABORT("FIXME\n");
}

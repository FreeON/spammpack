/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "chunk_block.h"
#include "chunk_tree.h"

#ifdef MIC_ALLOC
#include "chunk_mic.h"
#endif

#ifdef BLOCK_BLAS
#include "lapack_interface.h"
#endif

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
#define DEBUG(message, ...) printf("[%s:%d (%s) thread %d DEBUG] " message, \
    __FILE__, __LINE__, __func__, omp_get_thread_num(), ##__VA_ARGS__)
#else
#define DEBUG(message, ...) printf("[%s:%d (%s) DEBUG] " message, \
    __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#endif
#else
#define DEBUG(message, ...) /* stripped DEBUG statement. */
#endif

/** A convenience macro for printing some info level message. */
#ifdef _OPENMP
#define INFO(message, ...) printf("[%s:%d (%s) thread %d] " message, \
    __FILE__, __LINE__, __func__, omp_get_thread_num(), ##__VA_ARGS__)
#else
#define INFO(message, ...) printf("[%s:%d (%s)] " message, \
    __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#endif

/** A convenience macro for print a fatal error message and terminating the
 * code. */
#define ABORT(message, ...) printf("[%s:%d (%s) FATAL] " message, \
    __FILE__, __LINE__, __func__, ##__VA_ARGS__); exit(1)

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

  /** The complexity (operation count) of the last chunk operation. This
   * counter gets reset at the beginning of a multiply() or an add()
   * operation.
   */
  size_t complexity;

  /** The offset from data[0] to the array of sub-matrices, skipping the tree,
   * going straight to the sub-matrices attached to the leaf nodes. */
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

  /** The children nodes are linked using relative offsets so we can pass the
   * chunk safely through message marshalling in between chares. This also
   * allows us to checkpoint a chunk without having to worry about relinking
   * the nodes after reading the chunk back from disk.
   *
   * The offset is based on the absolute pointer to the node, i.e. the
   * children or the matrix is given relative to the current node. We can not
   * use a reference outside of the node struct since we wouldn't have access
   * to it inside the struct. And we really want to preserve the nice
   * recursive programming model of SpAMM insde the chunk.
   */
  union
  {
    /** The children nodes (cast as offsets based on pointer to this struct). */
    intptr_t child_offset[4];

    /** The submatrix for leaf nodes (cast as offset from the leaf node). */
    intptr_t matrix_offset;
  }
  offset;
};

/** A convenience function to get 2^i, where i is integer.
 *
 * @param i The integer.
 *
 * @return 2 raised to the power of i.
 */
static inline int
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
#ifdef MIC_ALLOC
  void *chunk = malloc_huge_pages(chunk_tree_sizeof(N_chunk, N_basic));
  memset(chunk, 0, chunk_tree_sizeof(N_chunk, N_basic));
#else
  void *chunk = calloc(chunk_tree_sizeof(N_chunk, N_basic), 1);
#endif
  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;

  ptr->chunksize = chunk_tree_sizeof(N_chunk, N_basic);
  ptr->i_lower = i_lower;
  ptr->j_lower = j_lower;
  ptr->N = N;
  ptr->N_chunk = N_chunk;
  ptr->N_basic = N_basic;
  ptr->depth = chunk_tree_get_depth(N_chunk, N_basic);
  ptr->complexity = 0;

  DEBUG("new tree chunk, N_chunk = %d, N_basic = %d, N = %d, depth = %d\n",
      N_chunk, N_basic, N, ptr->depth);

  DEBUG("allocating new chunk\n");
  DEBUG("ptr = %p\n", ptr);
  DEBUG("ptr->data = %p\n", ptr->data);
  DEBUG("ptr->data + chunksize = %p\n", (void*) ((intptr_t) ptr + ptr->chunksize));
  DEBUG("sizeof(struct chunk_tree_node_t) = 0x%lx\n", sizeof(struct chunk_tree_node_t));
  DEBUG("sizeof(submatrix) = 0x%lx\n", SQUARE(N_basic)*sizeof(double));

  /* Fill the chunk with a matrix tree. */
  DEBUG("linking tree nodes...\n");
  intptr_t tier_offset = 0;
  for(int tier = 0; tier < ptr->depth; tier++)
  {
    intptr_t next_tier_offset = tier_offset + ipow2(2*tier)*sizeof(struct chunk_tree_node_t);

    for(int i = 0; i < ipow2(tier); i++)
    {
      for(int j = 0; j < ipow2(tier); j++)
      {
        intptr_t node_offset = ROW_MAJOR(i, j, ipow2(tier))*sizeof(struct chunk_tree_node_t);
        struct chunk_tree_node_t *node = (struct chunk_tree_node_t*)
          ((intptr_t) ptr->data + tier_offset + node_offset);

        DEBUG("%d:node(%d,%d) at %p\n", tier, i, j, node);

        /* Initialize node. */
        node->N_basic = N_basic;

        for(int i_child = 0; i_child < 2; i_child++)
        {
          for(int j_child = 0; j_child < 2; j_child++)
          {
            /* Offset relative to node. */
            intptr_t child_offset = next_tier_offset
                + ROW_MAJOR(2*i+i_child, 2*j+j_child,
                    ipow2(tier+1)) * sizeof(struct chunk_tree_node_t)
                - (tier_offset + node_offset);

            /* The child_offset is relative to node. */
            node->offset.child_offset[ROW_MAJOR(i_child, j_child, 2)] = child_offset;
            DEBUG("%d:child(%d,%d) -> %d:node(%d,%d) at %p\n",
                tier, i_child, j_child, tier+1, 2*i+i_child, 2*j+j_child,
                (void*) ((intptr_t) node + child_offset));
          }
        }
      }
    }

    tier_offset = next_tier_offset;
  }

  /* The offset directly into the sub-matrix array of the leaf nodes. */
  intptr_t leaf_tier_size = ipow2(2*ptr->depth)*sizeof(struct chunk_tree_node_t);
  ptr->submatrix_offset = tier_offset + leaf_tier_size;
  DEBUG("submatrix_offset = %ld\n", ptr->submatrix_offset);

  /* Link leaf nodes and store pointers to the basic sub-matrices. */
  DEBUG("linking leaf nodes and submatrices...\n");
  for(int i = 0; i < ipow2(ptr->depth); i++)
  {
    for(int j = 0; j < ipow2(ptr->depth); j++)
    {
      /* Pointer to leaf node. */
      intptr_t node_offset = ROW_MAJOR(i, j, ipow2(ptr->depth))*sizeof(struct chunk_tree_node_t);
      struct chunk_tree_node_t *node = (struct chunk_tree_node_t*)
        ((intptr_t) ptr->data + tier_offset + node_offset);

      /* Initialize node. */
      node->N_basic = N_basic;

#ifdef _OPENMP
      omp_init_lock(&node->matrix_lock);
#endif

      size_t matrix_offset = ROW_MAJOR(i, j, ipow2(ptr->depth))*SQUARE(N_basic)*sizeof(double);
      node->offset.matrix_offset = leaf_tier_size + matrix_offset - node_offset;

      DEBUG("%d:node(%d,%d) at %p, submatrix at %p\n",
          ptr->depth, i, j, node, &(*(double*) ((intptr_t) node + node->offset.matrix_offset)));
    }
  }

  DEBUG("chunk ends at %p, done\n", (void *) ((intptr_t) ptr->data
        + ptr->submatrix_offset
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

  double *submatrix_ptr = (double*) ((intptr_t) chunk->data
      + chunk->submatrix_offset
      + ROW_MAJOR(i, j, ipow2(chunk->depth))*SQUARE(chunk->N_basic)*sizeof(double));

  DEBUG("submatrix(%d,%d) at %p\n", i, j, submatrix_ptr);

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
              node->norm_2 += SQUARE(submatrix[COLUMN_MAJOR(i_basic, j_basic, chunk->N_basic)]);
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
          submatrix_ptr[COLUMN_MAJOR(i_basic, j_basic, ptr->N_basic)] =
            A[COLUMN_MAJOR(i*ptr->N_basic+i_basic, j*ptr->N_basic+j_basic, ptr->N_chunk)];
        }
      }
    }
  }

#ifdef DEBUG_OUTPUT_MATRIX
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
          printf(" % 1.3f", submatrix_ptr[COLUMN_MAJOR(i_basic, j_basic, ptr->N_basic)]);
        }
        printf("\n");
      }
    }
  }
#endif

  /** Update norms. */
  chunk_tree_update_norm(chunk);
}

/** Return the matrix size of the chunk matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_chunk.
 */
int
chunk_tree_get_N_chunk (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_tree_t*) chunk)->N_chunk;
}

/** Return the matrix size of the basic matrix.
 *
 * @param chunk The chunk.
 *
 * @return The size N_basic.
 */
int
chunk_tree_get_N_basic (void *const chunk)
{
  assert(chunk != NULL);
  return ((struct chunk_tree_t*) chunk)->N_basic;
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

/** Get the square of the matrix norm of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The square of the Frobenius norm.
 */
double
chunk_tree_get_norm_2 (const void *const chunk)
{
  assert(chunk != NULL);
  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;
  struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;
  return root->norm_2;
}

/** Print a node.
 *
 * @param tier The tier.
 * @param depth The depth.
 * @param node The node.
 */
void
chunk_tree_print_node (const int tier,
    const int depth,
    const struct chunk_tree_node_t *node)
{
  assert(node != NULL);

  INFO("(%p) norm_2  = %e\n", node, node->norm_2);
  INFO("(%p) N_basic = %d\n", node, node->N_basic);

  if(tier == depth)
  {
    INFO("(%p) submatrix at %p\n", node, (void*) ((intptr_t) node + node->offset.matrix_offset));
    double *submatrix = (double*) ((intptr_t) node + node->offset.matrix_offset);
    for(int i = 0; i < node->N_basic; i++)
    {
      INFO("(%p) ", node);
      for(int j = 0; j < node->N_basic; j++)
      {
        printf(" % 1.3f", submatrix[COLUMN_MAJOR(i, j, node->N_basic)]);
      }
      printf("\n");
    }
  }

  else
  {
    INFO("(%p) children at [ %p %p\n",   node,
        (void*) ((intptr_t) node + node->offset.child_offset[0]),
        (void*) ((intptr_t) node + node->offset.child_offset[1]));
    INFO("(%p)               %p %p ]\n", node,
        (void*) ((intptr_t) node + node->offset.child_offset[2]),
        (void*) ((intptr_t) node + node->offset.child_offset[3]));

    for(int i = 0; i < 4; i++)
    {
      chunk_tree_print_node(tier+1, depth,
          (struct chunk_tree_node_t*) ((intptr_t) node + node->offset.child_offset[i]));
    }
  }
}

/** Print a chunk.
 *
 * @param chunk The chunk.
 * @param format The format string. This follows the printf() approach.
 */
void
chunk_tree_print (const void *const chunk,
    const char *const format, ...)
{
  INFO("Chunk -->\n");
  INFO("chunk at %p\n", chunk);

  if(chunk != NULL)
  {
    struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;

    INFO("chunksize        = %ld\n", ptr->chunksize);
    INFO("i_lower          = %d\n", ptr->i_lower);
    INFO("j_lower          = %d\n", ptr->j_lower);
    INFO("N                = %d\n", ptr->N);
    INFO("N_chunk          = %d\n", ptr->N_chunk);
    INFO("N_basic          = %d\n", ptr->N_basic);
    INFO("depth            = %d\n", ptr->depth);
    INFO("submatrix_offset = 0x%lx\n", ptr->submatrix_offset);

    struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;

    chunk_tree_print_node(0, ptr->depth, root);
  }
  INFO("<-- Chunk\n");
}

/** Add two nodes.
 *
 * @f[ A \leftarrow \alpha A + \beta B @f]
 *
 * @param tier The tier.
 * @param depth The depth.
 * @param alpha The scalar alpha.
 * @param A The node A.
 * @param beta The scalar beta.
 * @param B The node B.
 */
void
chunk_tree_add_node (const int tier,
    const int depth,
    const double alpha,
    struct chunk_tree_node_t *const A,
    const double beta,
    const struct chunk_tree_node_t *const B)
{
  assert(A != NULL);
  assert(B != NULL);

  if(tier == depth)
  {
    double *A_matrix = (double*) ((intptr_t) A + A->offset.matrix_offset);
    double *B_matrix = (double*) ((intptr_t) B + B->offset.matrix_offset);
    for(int i = 0; i < SQUARE(A->N_basic); i++)
    {
      A_matrix[i] = alpha*A_matrix[i] + beta*B_matrix[i];
    }
  }

  else
  {
    for(int i = 0; i < 4; i++)
    {
      struct chunk_tree_node_t *A_child = (struct chunk_tree_node_t*)
        ((intptr_t) A + A->offset.child_offset[i]);
      struct chunk_tree_node_t *B_child = (struct chunk_tree_node_t*)
        ((intptr_t) B + B->offset.child_offset[i]);

      chunk_tree_add_node(tier+1, depth, alpha, A_child, beta, B_child);
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
chunk_tree_add (const double alpha, void *const A,
    const double beta, const void *const B)
{
  assert(A != NULL);
  assert(B != NULL);

  struct chunk_tree_t *A_ptr = (struct chunk_tree_t*) A;
  struct chunk_tree_t *B_ptr = (struct chunk_tree_t*) B;

  struct chunk_tree_node_t *A_root = (struct chunk_tree_node_t*) A_ptr->data;
  struct chunk_tree_node_t *B_root = (struct chunk_tree_node_t*) B_ptr->data;

  DEBUG("adding two chunks\n");

#ifdef DEBUG_OUTPUT_MATRIX
  chunk_tree_print(A, "chunk A\n");
  chunk_tree_print(B, "chunk B\n");
#endif

  chunk_tree_add_node(0, A_ptr->depth, alpha, A_root, beta, B_root);
  chunk_tree_update_norm(A);
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
 * @param symbolic_only Only go through the symbolic part, i.e. don't multiply the
 * basic blocks. Used only for debugging to get a handle on the performance of
 * the tree.
 * @param complexity The complexity counter array.
 * @param product_index The linear index in the product space.
 */
void
chunk_tree_multiply_node (const double tolerance_2,
    const int tier,
    const int depth,
    const struct chunk_tree_node_t *const A,
    const struct chunk_tree_node_t *const B,
    struct chunk_tree_node_t *const C,
    const short symbolic_only,
    size_t *const complexity,
    const int product_index)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  DEBUG("%d(%d) A at %p, B at %p, C at %p\n", tier, depth, A, B, C);

  DEBUG("symbolic_only = %d\n", symbolic_only);

  if(tier == depth)
  {
    DEBUG("%d(%d) A at %p, B at %p, C at %p\n", tier, depth, A, B, C);

    if(!symbolic_only)
    {
      double *A_submatrix = (double*) ((intptr_t) A + A->offset.matrix_offset);
      double *B_submatrix = (double*) ((intptr_t) B + B->offset.matrix_offset);
      double *C_submatrix = (double*) ((intptr_t) C + C->offset.matrix_offset);

      DEBUG("multiplying submatrix, A at %p, B at %p, C at %p\n",
          &(*A_submatrix), &(*B_submatrix), &(*C_submatrix));

#ifdef MEASURE_COMPLEXITY
      complexity[product_index]++;
#endif

      //INFO("A[0][0] = % 1.3f\n", A_submatrix[0]);
      //INFO("B[0][0] = % 1.3f\n", B_submatrix[0]);
      //INFO("C[0][0] = % 1.3f\n", C_submatrix[0]);

#ifdef _OPENMP
      omp_set_lock(&C->matrix_lock);
#endif

#if defined(BLOCK_MULTIPLY)
      chunk_block_multiply(A_submatrix, B_submatrix, C_submatrix, A->N_basic);
#elif defined(BLOCK_BLAS)
      {
        double alpha = 1.0;
        double beta = 1.0;
        DGEMM("N", "N", &A->N_basic, &A->N_basic, &A->N_basic, &alpha,
            A_submatrix, &A->N_basic, B_submatrix, &A->N_basic, &beta,
            C_submatrix, &A->N_basic);
      }
#else
#error unknown block multiply implementation.
#endif

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
#pragma omp task default(none) firstprivate(i, j, k) untied if(tier <= CHUNK_TREE_MAX_TIER)
          {
            DEBUG("(%d,%d,%d) starting new task\n", i, j, k);

            struct chunk_tree_node_t *A_child = (struct chunk_tree_node_t*)
              ((intptr_t) A + A->offset.child_offset[ROW_MAJOR(i, k, 2)]);
            struct chunk_tree_node_t *B_child = (struct chunk_tree_node_t*)
              ((intptr_t) B + B->offset.child_offset[ROW_MAJOR(k, j, 2)]);
            struct chunk_tree_node_t *C_child = (struct chunk_tree_node_t*)
              ((intptr_t) C + C->offset.child_offset[ROW_MAJOR(i, j, 2)]);

            DEBUG("(%d,%d,%d) A at %p, B at %p, C at %p\n", i, j, k, A_child, B_child, C_child);

            if(A_child->norm_2*B_child->norm_2 > tolerance_2)
            {
#ifdef MEASURE_COMPLEXITY
              int new_product_index = (product_index << 3) | (i << 2) | (j << 1) | k;
#else
              int new_product_index = product_index;
#endif
              chunk_tree_multiply_node(tolerance_2, tier+1, depth, A_child,
                  B_child, C_child, symbolic_only, complexity, new_product_index);
            }

            else
            {
              DEBUG("skipping product, A_norm*B_norm = %e\n", A_child->norm_2*B_child->norm_2);
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
 * @param symbolic_only Only go through the symbolic part, i.e. don't multiply the
 * basic blocks. Used only for debugging to get a handle on the performance of
 * the tree.
 */
void
chunk_tree_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C,
    const short symbolic_only)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  struct chunk_tree_t *A_ptr = (struct chunk_tree_t*) A;
  struct chunk_tree_t *B_ptr = (struct chunk_tree_t*) B;
  struct chunk_tree_t *C_ptr = (struct chunk_tree_t*) C;

  /* Reset C. */
  chunk_tree_set_zero(C_ptr);
  C_ptr->complexity = 0;

  double tolerance_2 = SQUARE(tolerance);

  struct chunk_tree_node_t *A_root = (struct chunk_tree_node_t*) A_ptr->data;
  struct chunk_tree_node_t *B_root = (struct chunk_tree_node_t*) B_ptr->data;
  struct chunk_tree_node_t *C_root = (struct chunk_tree_node_t*) C_ptr->data;

  DEBUG("symbolic_only = %d\n", symbolic_only);

  DEBUG("%dx%d blocked matrix, potentially %d products to consider, max tier = %d\n",
      ipow2(A_ptr->depth), ipow2(A_ptr->depth), CUBE(ipow2(A_ptr->depth), CHUNK_TREE_MAX_TIER));

  DEBUG("SpAMM tolerance = %e\n", tolerance);

#ifdef MEASURE_COMPLEXITY
  size_t *complexity = calloc(CUBE(A_ptr->N_chunk/A_ptr->N_basic), sizeof(size_t));
#else
  void *complexity = NULL;
#endif

  if(A_root->norm_2*B_root->norm_2 > tolerance_2)
  {
#pragma omp parallel
    {
#pragma omp master
      {
#ifdef _OPENMP
        DEBUG("running on %d OpenMP threads\n", omp_get_num_threads());
#else
        DEBUG("running in serial\n");
#endif

#pragma omp task untied
        {
          chunk_tree_multiply_node(tolerance_2, 0, A_ptr->depth, A_root,
              B_root, C_root, symbolic_only, complexity, 0);
        }
#pragma omp taskwait
      }
    }
  }

#ifdef MEASURE_COMPLEXITY
  size_t product_complexity = 0;
#pragma omp parallel for reduction(+:product_complexity)
  for(int i = 0; i < CUBE(A_ptr->N_chunk/A_ptr->N_basic); i++)
  {
    product_complexity += complexity[i];
  }
  free(complexity);

  C_ptr->complexity += product_complexity;
  DEBUG("product complexity = %ld out of %d, complexity ratio = %1.3f\n",
      product_complexity, CUBE(A_ptr->N_chunk/A_ptr->N_basic),
      product_complexity/(double) CUBE(A_ptr->N_chunk/A_ptr->N_basic));
#endif
}

/** Convert a chunk to a dense matrix.
 *
 * @param chunk The chunk.
 *
 * @return The dense matrix in column-major order. This matrix is sized to
 * blocksize, which is stored in the chunk.
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
            submatrix[COLUMN_MAJOR(i_basic, j_basic, ptr->N_basic)];
        }
      }
    }
  }

  return A;
}

/** Scale a node, i.e. multiply it by a scalar.
 *
 * @f[ A \leftarrow \alpha A @f]
 *
 * @param alpha The scalar alpha.
 * @param tier The tier.
 * @param depth The depth.
 * @param node The node.
 */
void
chunk_tree_scale_node (const double alpha,
    const int tier,
    const int depth,
    struct chunk_tree_node_t *node)
{
  assert(node != NULL);

  DEBUG("scaling node with alpha = %e\n", alpha);

  node->norm_2 *= SQUARE(alpha);

  if(tier == depth)
  {
    double *matrix = (double*) ((intptr_t) node + node->offset.matrix_offset);
    for(int i = 0; i < SQUARE(node->N_basic); i++)
    {
      matrix[i] *= alpha;
    }
  }

  else
  {
    for(int i = 0; i < 4; i++)
    {
      struct chunk_tree_node_t *child = (struct chunk_tree_node_t*)
        ((intptr_t) node + node->offset.child_offset[i]);
      chunk_tree_scale_node(alpha, tier+1, depth, child);
    }
  }
}

/** Scale a chunk, i.e. multiply it by a scalar.
 *
 * @f[ A \leftarrow \alpha A @f]
 *
 * @param alpha The scalar alpha.
 * @param chunk The chunk.
 */
void
chunk_tree_scale (const double alpha, void *const chunk)
{
  assert(chunk != NULL);

  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;
  struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;

  DEBUG("scaling tree with alpha = %e\n", alpha);
  chunk_tree_scale_node(alpha, 0, ptr->depth, root);
}

/** Add a scaled identity matrix to a chunk.
 *
 * @f[ A \leftarrow A + \alpha \times \mathrm{Id} @f]
 *
 * @param alpha The scalar alpha.
 * @param node The node.
 */
void
chunk_tree_add_identity_node (const int tier,
    const int depth,
    const double alpha,
    struct chunk_tree_node_t *const node)
{
  assert(node != NULL);

  if(tier == depth)
  {
    double *matrix = (double*) ((intptr_t) node + node->offset.matrix_offset);
    for(int i = 0; i < node->N_basic; i++)
    {
      matrix[COLUMN_MAJOR(i, i, node->N_basic)] += alpha;
    }
  }

  else
  {
    for(int i = 0; i < 2; i++)
    {
      struct chunk_tree_node_t *child = (struct chunk_tree_node_t*)
        ((intptr_t) node + node->offset.child_offset[ROW_MAJOR(i, i, 2)]);
      chunk_tree_add_identity_node(tier+1, depth, alpha, child);
    }
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
chunk_tree_add_identity (const double alpha, void *const chunk)
{
  assert(chunk != NULL);

  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;
  struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;

  chunk_tree_add_identity_node(0, ptr->depth, alpha, root);
  chunk_tree_update_norm(chunk);
}

/** Get the trace of a node.
 *
 * @param node The node.
 *
 * @return The trace of the node.
 */
double
chunk_tree_trace_node (const int tier,
    const int depth,
    const struct chunk_tree_node_t *const node)
{
  assert(node != NULL);

  double trace = 0;

  if(tier == depth)
  {
    double *matrix = (double*) ((intptr_t) node + node->offset.matrix_offset);
    for(int i = 0; i < node->N_basic; i++)
    {
      trace += matrix[COLUMN_MAJOR(i, i, node->N_basic)];
    }
  }

  else
  {
    for(int i = 0; i < 2; i++)
    {
      struct chunk_tree_node_t *child = (struct chunk_tree_node_t*)
        ((intptr_t) node + node->offset.child_offset[ROW_MAJOR(i, i, 2)]);
      trace += chunk_tree_trace_node(tier+1, depth, child);
    }
  }

  return trace;
}

/** Get the trace of a chunk.
 *
 * @param chunk The chunk.
 *
 * @return The trace of the chunk.
 */
double
chunk_tree_trace (const void *const chunk)
{
  assert(chunk != NULL);

  struct chunk_tree_t *ptr = (struct chunk_tree_t*) chunk;
  struct chunk_tree_node_t *root = (struct chunk_tree_node_t*) ptr->data;

  return chunk_tree_trace_node(0, ptr->depth, root);
}

/** Delete the chunk.
 *
 * @param chunk The chunk.
 */
void
chunk_tree_delete (void **const chunk)
{
  /* Ignore OpenMP locks for now. */
#ifdef MIC_ALLOC
  free_huge_pages(*chunk);
#else
  free(*chunk);
#endif
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
chunk_tree_get_complexity (const void *const chunk)
{
  const struct chunk_tree_t *ptr = chunk;
  return ptr->complexity;
}

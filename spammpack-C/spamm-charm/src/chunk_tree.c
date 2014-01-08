/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

  INFO("ptr = %p\n", ptr);
  INFO("ptr->data = %p\n", ptr->data);
  INFO("ptr->data + chunksize = %p\n", (void*) ((intptr_t) ptr + ptr->chunksize));
  INFO("sizeof(struct chunk_tree_node_t) = 0x%lx\n", sizeof(struct chunk_tree_node_t));
  INFO("sizeof(submatrix) = 0x%lx\n", SQUARE(N_basic)*sizeof(double));

  /* Fill the chunk with a matrix tree. */
  INFO("linking tree nodes...\n");
  intptr_t tier_ptr = (intptr_t) ptr->data;
  for(int tier = 0; tier < ptr->depth; tier++)
  {
    intptr_t next_tier_ptr = tier_ptr + ipow2(2*tier)*sizeof(struct chunk_tree_node_t);

    for(int i = 0; i < ipow2(tier); i++)
    {
      for(int j = 0; j < ipow2(tier); j++)
      {
        size_t offset = (i*ipow2(tier)+j)*sizeof(struct chunk_tree_node_t);
        struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) (tier_ptr + offset);

        INFO("%d:node(%d,%d) at %p\n", tier, i, j, node);

        for(int i_child = 0; i_child < 2; i_child++)
        {
          for(int j_child = 0; j_child < 2; j_child++)
          {
            size_t child_offset = ((2*i+i_child)*ipow2(tier+1)+2*j+j_child)*sizeof(struct chunk_tree_node_t);
            struct chunk_tree_node_t *child = (struct chunk_tree_node_t*) (next_tier_ptr + child_offset);
            node->data.child[i_child*2+j_child] = child;
            INFO("%d:child(%d,%d) -> %d:node(%d,%d) at %p\n",
                tier, i_child, j_child, tier+1, 2*i+i_child, 2*j+j_child,
                child);
          }
        }
      }
    }

    tier_ptr = next_tier_ptr;
  }

  /* Store pointers to the basic sub-matrices. */
  INFO("linking submatrices...\n");
  intptr_t submatrix_ptr = tier_ptr + ipow2(2*ptr->depth)*sizeof(struct chunk_tree_node_t);
  for(int i = 0; i < ipow2(ptr->depth); i++)
  {
    for(int j = 0; j < ipow2(ptr->depth); j++)
    {
      size_t offset = (i*ipow2(ptr->depth)+j)*sizeof(struct chunk_tree_node_t);
      size_t matrix_offset = (i*ipow2(ptr->depth)+j)*SQUARE(N_basic)*sizeof(double);
      struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) (tier_ptr + offset);
      node->data.matrix = (double*) (submatrix_ptr + matrix_offset);
      INFO("%d:node(%d,%d) at %p, submatrix at %p\n", ptr->depth, i, j, node, &(*node->data.matrix));
    }
  }

  INFO("chunk ends at %p, done\n", (void *) (submatrix_ptr
        + ipow2(2*ptr->depth)*SQUARE(N_basic)*sizeof(double)));

  return chunk;
}

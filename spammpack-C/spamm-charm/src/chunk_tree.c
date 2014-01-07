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

  else if(i == 1)
  {
    return 1 << i;
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
  if(ipow2(depth) != N_chunk)
  {
    printf("logic failure\n");
    exit(1);
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
  size_t matrix_size = ipow2(2*chunk_tree_get_depth(N_chunk, N_basic))*N_basic*N_basic*sizeof(double);
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

  /* Fill the chunk with a matrix tree. */
  for(int tier = 0; tier < ptr->depth; tier++)
  {
    for(int i = 0; i < ipow2(tier); i++)
    {
      for(int j = 0; j < ipow2(tier); j++)
      {
        struct chunk_tree_node_t *node = (struct chunk_tree_node_t*) ((intptr_t) ptr->data
            + (i*ipow2(tier)+j)*ipow2(2*tier)*sizeof(struct chunk_tree_node_t));

        for(int i_child = 0; i_child < 2; i_child++)
        {
          for(int j_child = 0; j_child < 2; j_child++)
          {
            struct chunk_tree_node_t *child = (struct chunk_tree_node_t*) ((intptr_t) ptr->data
                + ((2*i+i_child)*ipow2(tier)+2*j+j_child)*ipow2(2*(tier+1))*sizeof(struct chunk_tree_node_t));
            node->data.child[i_child*2+j_child] = child;
          }
        }
      }
    }
  }

  /* Store pointers to the basic sub-matrices. */

  return chunk;
}

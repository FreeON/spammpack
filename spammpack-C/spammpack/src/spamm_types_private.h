/** @file */

#ifndef __SPAMM_TYPES_PRIVATE_H
#define __SPAMM_TYPES_PRIVATE_H

#include "config.h"

/** The matrix. */
struct spamm_matrix_t
{
  /** The number of dimensions of this matrix (a 1-dimensional matrix is a
   * vector, a 2-dimensional matrix is a matrix, and so on...).
   */
  unsigned int number_dimensions;

  /** Number of rows/columns in this matrix. The length of this array is given
   * by number_dimensions. */
  unsigned int *N;

  /** Number of rows and columns in the padded matrix. */
  unsigned int N_padded;

  /** The block size at which the SpAMM condition is applied. */
  unsigned int N_block;

  /** The total depth of the tree. */
  unsigned int depth;

  /** The tier of the contiguous matrix chunks. At this tier the submatrix is
   * stored in SpAMM chunks, and accessed hierarchically or through a linear
   * subtree. */
  unsigned int contiguous_tier;

  /** Use a linear submatrix representation or a hierarchical one. A value of
   * 0 implies hierarchical, while anything else (but we use 1) implies
   * linear. */
  short use_linear_tree;

  /** The kernel tier. At this tier the stream kernel processes matrices. */
  unsigned int kernel_tier;

  union tree_t
  {
    /** The root node of the recursive tree. */
    struct spamm_recursive_node_t *recursive_tree;

    /** The hashed tree (if we deal with a hybrid matrix tree). */
    struct spamm_hashed_t *hashed_tree;

    /** The SpAMM chunk. */
    spamm_chunk_t *chunk;
  }
  tree;
};

/** The hashed matrix type. */
struct spamm_hashed_t
{
  /** The layout of the basic matrix blocks on the kernel tier. */
  enum spamm_layout_t layout;

  /** The lower value of the row index. */
  unsigned int M_lower;

  /** The upper value of the row index. */
  unsigned int M_upper;

  /** The lower value of the column index. */
  unsigned int N_lower;

  /** The upper value of the column index. */
  unsigned int N_upper;

  /** The tier this linear tree is attached to. */
  unsigned int tier;

  /** Tree depth. */
  unsigned int depth;

  /** The kernel tier. */
  unsigned int kernel_tier;

  /** The hashtables for access to each tier. */
  struct spamm_hashtable_t **tier_hashtable;
};

/** A node in the matrix tree. */
struct spamm_hashed_node_t
{
  /** The tier. */
  unsigned int tier;

  /** The linear index of this node in 2D, i.e. in matrix space. */
  unsigned int index_2D;

  /** The norm of this block. */
  float norm;

  /** The square of the norm of this block. */
  float norm2;
};

/** A node at the kernel tier. */
struct spamm_hashed_data_t
{
  /** The layout of the basic matrix blocks on the kernel tier. */
  enum spamm_layout_t layout;

  /** The tier. */
  unsigned int tier;

  /** The linear index of this node in 2D, i.e. in matrix space. */
  unsigned int index_2D;

  /** The norm of this matrix block. */
  float node_norm;

  /** The square of the norm of this matrix block. */
  float node_norm2;

  /** The norms of the basic block matrices.
   *
   * The norms are the Frobenius norms for the basic #SPAMM_N_BLOCK x
   * #SPAMM_N_BLOCK matrix blocks. Since they are #SPAMM_N_KERNEL_BLOCKED x
   * #SPAMM_N_KERNEL_BLOCKED of those, the norms are stored in row-major order
   * in a #SPAMM_N_KERNEL_BLOCKED x #SPAMM_N_KERNEL_BLOCKED matrix. */
  float __attribute__ ((aligned (8))) norm[SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED];

  /** The square of the norms of the basic block matrices. */
  float __attribute__ ((aligned (8))) norm2[SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED];

  /** The matrix data.
   *
   * The matrix elements are arranged on a row-major ordered
   * #SPAMM_N_KERNEL_BLOCKED x #SPAMM_N_KERNEL_BLOCKED grid of basic *
   * #SPAMM_N_BLOCK x #SPAMM_N_BLOCK matrix blocks. The total size of
   * block_dense is therefore (#SPAMM_N_KERNEL_BLOCKED * #SPAMM_N_BLOCK) x
   * (#SPAMM_N_KERNEL_BLOCKED * #SPAMM_N_BLOCK). */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense[SPAMM_N_KERNEL*SPAMM_N_KERNEL];

  /** A store buffer. */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense_store[SPAMM_N_KERNEL*SPAMM_N_KERNEL];

  /** The transpose of block_dense.  */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense_transpose[SPAMM_N_KERNEL*SPAMM_N_KERNEL];

  /** The matrix data (dilated by 4 for SSE). */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense_dilated[SPAMM_N_KERNEL*SPAMM_N_KERNEL*4];
};

/** A recursive node. */
struct
spamm_recursive_node_t
{
  /** The norm of this block. */
  float norm;

  /** The square of the norm of this block. */
  float norm2;

  union children
  {
    /** The children nodes. */
    struct spamm_recursive_node_t **child;

    /** The hashed tree (if we deal with a hybrid matrix tree). */
    struct spamm_hashed_t *hashed_tree;

    /** The matrix data (if this is a leaf node) of size N_contiguous x
     * N_contiguous. */
    spamm_chunk_t *chunk;
  }
  tree;
};

/** The basic information in a stream element. */
struct spamm_multiply_stream_t
{
  /** A pointer to the kernel tier matrix node of A. */
  struct spamm_hashed_data_t *A;

  /** A pointer to the kernel tier matrix node of B. */
  struct spamm_hashed_data_t *B;

  /** A pointer to the kernel tier matrix node of C. */
  struct spamm_hashed_data_t *C;
};

#endif

/** @file */

#ifndef __SPAMM_TYPES_PRIVATE_H
#define __SPAMM_TYPES_PRIVATE_H

#include "spamm_config.h"

/** The matrix.
 */
struct spamm_matrix_t
{
  /** Number or rows in this matrix. */
  unsigned int M;

  /** Number of columns in this matrix. */
  unsigned int N;

  /** Number of rows and columns in the padded matrix. */
  unsigned int N_padded;

  /** The total depth of the tree. */
  unsigned int depth;

  /** The number of tiers from the bottom that are stored in the hashed
   * format. */
  unsigned int number_hashed_tiers;

  /** The kernel tier, i.e. the tier at which the dense kernel matrix blocks
   * will be processed. */
  unsigned int kernel_tier;

  /** The blocksize, i.e. the matrix size at the kernel tier. */
  unsigned int blocksize;

  /** The root node. */
  struct spamm_recursive_node_t *root;

  /** The hashed tree (if we deal with a hybrid matrix tree). */
  struct spamm_hashed_t *hashed_tree;
};

/** The hashed matrix type.
 */
struct spamm_hashed_t
{
  /** The layout of the basic matrix blocks on the kernel tier. */
  enum spamm_layout_t layout;

  /** Number or rows in this matrix. */
  unsigned int M;

  /** Number of columns in this matrix. */
  unsigned int N;

  /** Number of rows and columns in the padded matrix. */
  unsigned int N_padded;

  /** Tree depth. */
  unsigned int depth;

  /** The kernel tier. */
  unsigned int kernel_tier;

  /** The hashtables for access to each tier. */
  struct spamm_hashtable_t **tier_hashtable;
};

/** A node in the matrix tree.
 */
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
struct spamm_data_t
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

  /** The upper tier norms. */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) norm_upper[2*2*2];

  /** The upper tier norms for the transpose. */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) norm_upper_transpose[2*2*2];

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

/** A recursive recursive SpAMM tree. */
struct
spamm_recursive_t
{
  /** Number or rows in this matrix. */
  unsigned int M;

  /** Number of columns in this matrix. */
  unsigned int N;

  /** Number of rows and columns in the padded matrix. */
  unsigned int N_padded;

  /** Tree depth. */
  unsigned int depth;

  /** The blocksize. */
  unsigned int blocksize;

  /** The root node. */
  struct spamm_recursive_node_t *root;
};

struct
spamm_recursive_node_t
{
  /** The lower value of the row index. */
  unsigned int M_lower;

  /** The upper value of the row index. */
  unsigned int M_upper;

  /** The lower value of the column index. */
  unsigned int N_lower;

  /** The upper value of the column index. */
  unsigned int N_upper;

  /** The tier. */
  unsigned int tier;

  /** The norm of this block. */
  float norm;

  /** The square of the norm of this block. */
  float norm2;

  /** The children nodes. */
  struct spamm_recursive_node_t *child[4];

  /** The blocksize. */
  unsigned int blocksize;

  /** The matrix data (if this is a leaf node). */
  float *data;

  /** The hashed tree (if we deal with a hybrid matrix tree). */
  struct spamm_hashed_t *hashed_tree;
};

/** The basic information in a stream element.
 */
struct spamm_multiply_stream_t
{
  /** A pointer to the kernel tier matrix node of A. */
  struct spamm_data_t *A;

  /** A pointer to the kernel tier matrix node of B. */
  struct spamm_data_t *B;

  /** A pointer to the kernel tier matrix node of C. */
  struct spamm_data_t *C;
};

#endif

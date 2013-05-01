/** @file */

#ifndef __SPAMM_TYPES_PRIVATE_H
#define __SPAMM_TYPES_PRIVATE_H

#include "config.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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

  /** The total depth of the tree. */
  unsigned int depth;

  /** The tier of the contiguous matrix chunks. At this tier the submatrix is
   * stored in SpAMM chunks. If use_linear_tree, then the spamm_hashed_* code
   * takes over and processes the chunks, if not, then the chunk is
   * interpreted has containing a submatrix of dimension N_contiguous x
   * N_contiguous. */
  unsigned int chunk_tier;

  /** Use a linear submatrix representation or a hierarchical one. A value of
   * 0 implies hierarchical, while anything else (but we use 1) implies
   * linear. */
  short use_linear_tree;

  /** The root node of the recursive tree. */
  struct spamm_recursive_node_t *recursive_tree;
};

/** A recursive node. */
struct
spamm_recursive_node_t
{
  /** The reference count, i.e. a count of how many non-zero nodes are
   * underneath this node.
   */
  int refcount;

  /** The norm of this block. */
  double norm;

  /** The square of the norm of this block. */
  double norm2;

#ifdef _OPENMP
  /** A lock. */
  omp_lock_t lock;
#endif

  union children
  {
    /** The children nodes. */
    struct spamm_recursive_node_t **child;

    /** The matrix data (if this is a leaf node) of size N_contiguous x
     * N_contiguous. */
    spamm_chunk_t *chunk;
  }
  tree;
};

#endif

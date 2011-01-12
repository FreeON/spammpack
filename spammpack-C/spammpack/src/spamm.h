/** @file */

#ifndef __SPAMM_H
#define __SPAMM_H

/* Include some configuration options. */
#include "spamm_config.h"
#include "spamm_error.h"
#include "spamm_hashtable.h"
#include "spamm_list.h"
#include "spamm_kernel.h"
#include "spamm_timer.h"
#include "spamm_convert.h"

/** The matrix type.
 */
struct spamm_t
{
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
struct spamm_node_t
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
   * #SPAMM_N_BLOCK matrix blocks. Since they are #SPAMM_N_KERNEL_BLOCK x
   * #SPAMM_N_KERNEL_BLOCK of those, the norms are stored in row-major order
   * in a #SPAMM_N_KERNEL_BLOCK x #SPAMM_N_KERNEL_BLOCK matrix. */
  float norm[SPAMM_N_KERNEL_BLOCK*SPAMM_N_KERNEL_BLOCK];

  /** The square of the norms of the basic block matrices. */
  float norm2[SPAMM_N_KERNEL_BLOCK*SPAMM_N_KERNEL_BLOCK];

  /** The matrix data.
   *
   * The matrix elements are arranged on a row-major ordered
   * #SPAMM_N_KERNEL_BLOCK x #SPAMM_N_KERNEL_BLOCK grid of basic *
   * #SPAMM_N_BLOCK x #SPAMM_N_BLOCK matrix blocks. The total size of
   * block_dense is therefore (#SPAMM_N_KERNEL_BLOCK * #SPAMM_N_BLOCK) x
   * (#SPAMM_N_KERNEL_BLOCK * #SPAMM_N_BLOCK). */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense[SPAMM_N_KERNEL*SPAMM_N_KERNEL];

  /** The matrix data (dilated by 4 for SSE). */
  float __attribute__ ((aligned (SPAMM_ALIGNMENT))) block_dense_dilated[SPAMM_N_KERNEL*SPAMM_N_KERNEL*4];
};

/* Function declarations. */

void *
spamm_allocate (size_t size);

int
spamm_check (const struct spamm_t *A);

void
spamm_delete (struct spamm_t **A);

void
spamm_delete_block (struct spamm_data_t **data);

void
spamm_delete_node (struct spamm_node_t **node);

unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j);

unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_column_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j);

//inline unsigned int
//spamm_dense_index (const unsigned int i, const unsigned int j)
//{
//  if (i >= SPAMM_N_KERNEL || i >= SPAMM_N_KERNEL)
//  {
//    fprintf(stderr, "i or j out of bounds\n");
//    exit(1);
//  }
//
//  return i*SPAMM_N_KERNEL+j;
//}

float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A);

unsigned int
spamm_index_2D (const unsigned int i, const unsigned int j);

void
spamm_index_2D_to_ij (const unsigned int index, unsigned int *i, unsigned int *j);

unsigned int
spamm_index_3D_0kj (const unsigned int k, const unsigned int j);

unsigned int
spamm_index_3D_ik0 (const unsigned int i, const unsigned int k);

unsigned int
spamm_index_3D_i0j_to_2D (const unsigned int index_3D_i0j);

unsigned int
spamm_index_3D_ikj_to_k (const unsigned int index_3D_ikj);

void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C);

unsigned int
spamm_number_nonzero (const struct spamm_t *A);

void
spamm_print_dense (const unsigned int M, const unsigned int N, const float *A);

struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N);

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index_2D);

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D);

void
spamm_print (const struct spamm_t *A);

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A);

void
spamm_uint_to_bin_string (const unsigned int i, char *result);

char *
spamm_version ();

#endif

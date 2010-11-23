#ifndef __SPAMM_H__
#define __SPAMM_H__

/* Include some configuration options. */
#include "spamm_config.h"
#include "spamm_kernel.h"

/* Include more header files. */
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

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
  GHashTable *tier_hashtable;
};

/** A node in the matrix tree.
 */
struct spamm_node_t
{
  /** The tier. */
  unsigned int tier;

  /** The linear index of this node in 2D, i.e. in matrix space. */
  unsigned int index_2D;

  /** The linear index of this node in 3D, i.e. in product space. */
  unsigned int index_3D_ik0;

  /** The linear index of this node in 3D, i.e. in product space. */
  unsigned int index_3D_0kj;
};

/** A node at the kernel tier. */
struct spamm_data_t
{
  /** The tier. */
  unsigned int tier;

  /** The linear index of this node in 2D, i.e. in matrix space. */
  unsigned int index_2D;

  /** The linear index of this node in 3D, i.e. in product space. */
  unsigned int index_3D_ik0;

  /** The linear index of this node in 3D, i.e. in product space. */
  unsigned int index_3D_0kj;

  /** The norms of the basic block matrices. */
  float norm[SPAMM_N_KERNEL_BLOCK*SPAMM_N_KERNEL_BLOCK];

  /** The matrix data. */
  float block_dense[SPAMM_N_KERNEL*SPAMM_N_KERNEL];
};

/* Function declarations. */

void
spamm_delete (struct spamm_t **A);

void
spamm_delete_block (struct spamm_data_t **data);

void
spamm_delete_node (struct spamm_node_t **node);

unsigned int
spamm_dense_index (const unsigned int i, const unsigned int j);

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

gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b);

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

void
spamm_multiply (const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C);

struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N);

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index_2D,
    const unsigned int index_3D_ik0, const unsigned int index_3D_0kj);

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D,
    const unsigned int index_3D_ik0, const unsigned int index_3D_0kj);

void
spamm_print (const struct spamm_t *A);

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A);

void
spamm_uint_to_bin_string (const unsigned int i, char *result);

#endif

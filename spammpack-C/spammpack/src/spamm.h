#ifndef __SPAMM_H__
#define __SPAMM_H__

/* Include some configuration options. */
#include "spamm_config.h"

/* Include more header files. */
#include <glib.h>

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
  GHashTable *tier;
};

/** A node in the matrix tree.
 */
struct spamm_node_t
{
  /** The tier. */
  unsigned int tier;

  /** The linear index of this node. */
  unsigned int index;
};

/** A node at the kernel tier. */
struct spamm_data_t
{
  /** The tier. */
  unsigned int tier;

  /** The linear index of this node. */
  unsigned int index;

  /** The matrix data. */
  float block_dense[SPAMM_N_KERNEL][SPAMM_N_KERNEL];
};

/* Function declarations. */

void
spamm_delete (struct spamm_t **A);

void
spamm_delete_node (struct spamm_node_t **node);

float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A);

gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b);

unsigned int
spamm_index_2D (const unsigned int i, const unsigned int j);

struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N);

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index);

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index);

void
spamm_print (const struct spamm_t *A);

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A);

#endif

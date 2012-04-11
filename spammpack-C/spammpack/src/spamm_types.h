/** @file */

#ifndef __SPAMM_TYPES_H
#define __SPAMM_TYPES_H

#include "spamm_config.h"

/** The layout type for the layout of the basic matrix blocks on the kernel
 * tier.
 */
enum spamm_layout_t
{
  /** Layout in row-major order. */
  row_major,

  /** Layout in column-major order. */
  column_major,

  /** Layout in Z-curve order. */
  Z_curve,

  /** Layout as a dense kernel block, in row-major order. */
  dense_column_major
};

/** The matrix type.
 */
struct spamm_t;

/** A node in the matrix tree.
 */
struct spamm_node_t;

/** A node at the kernel tier. */
struct spamm_data_t;

/** A recursive recursive SpAMM tree. */
struct
spamm_recursive_t;

struct
spamm_recursive_node_t;

#endif

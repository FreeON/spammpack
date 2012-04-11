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

struct spamm_multiply_stream_t;
struct spamm_t;
struct spamm_node_t;
struct spamm_data_t;
struct spamm_recursive_t;
struct spamm_recursive_node_t;

#endif

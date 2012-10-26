#include "config.h"
#include "spamm.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define BITFIELD_SIZE 15

/** Convert an n-dimensions index tuple to a linear index.
 *
 * @param number_dimensions The number of dimensions.
 * @param i The array of indices.
 *
 * @return The linear index.
 */
unsigned int
spamm_index_linear (const unsigned int number_dimensions,
    const unsigned int *const i)
{
  short bit_index;
  int dim;
  unsigned int getmask = 1;
  unsigned int setmask = 1;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < sizeof(unsigned int)*8; bit_index += number_dimensions)
  {
    for (dim = number_dimensions-1; dim >= 0 && bit_index+number_dimensions-1-dim < sizeof(unsigned int)*8; dim--)
    {
      if (i[dim] & getmask)
      {
        index |= setmask;
      }
      setmask <<= 1;
    }
    getmask <<= 1;
  }

  return index;
}

/** Convert index pair (i,j) into linear 2D index.
 *
 * @param i The row index.
 * @param j The column index.
 *
 * @return The linear 2D matrix index.
 */
unsigned int
spamm_index_2D (const unsigned int i, const unsigned int j)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 1;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if ((j & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;

    if ((i & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;
    getmask <<= 1;
  }

  return index;
}

/** Convert a linear 2D index to an index pair (i,j). If any of the two
 * indices i and j point to NULL, then this index is ignored and not returned.
 *
 * @param index The 2D linear index.
 * @param i The row index.
 * @param j The column index.
 */
void
spamm_index_2D_to_ij (const unsigned int index, unsigned int *i, unsigned int *j)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 1;

  if (i != NULL)
  {
    *i = 0;
  }

  if (j != NULL)
  {
    *j = 0;
  }

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if (j != NULL)
    {
      if ((index & getmask) != 0)
      {
        *j |= setmask;
      }
    }
    getmask <<= 1;

    if (i != NULL)
    {
      if ((index & getmask) != 0)
      {
        *i |= setmask;
      }
    }
    setmask <<= 1;
    getmask <<= 1;
  }
}

/** Convert an index pair into a convolution space 3D linear index. The index
 * pair is interpreted as \f$A(k,j)\f$, i.e. this function is meant to be
 * applied to matrix \f$B\f$ in the product \f$A \times B\f$.
 *
 * @param k The row index.
 * @param j The column index.
 *
 * @return The 3D linear index with 0kj ordering.
 */
unsigned int
spamm_index_3D_0kj (const unsigned int k, const unsigned int j)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 1;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if ((j & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;

    if ((k & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 2;
    getmask <<= 1;
  }

  return index;
}

/** Convert an index pair into a convolution space 3D linear index. The index
 * pair is interpreted as \f$A(i,k)\f$, i.e. this function is meant to be
 * applied to matrix \f$A\f$ in the product \f$A \times B\f$.
 *
 * @param i The row index.
 * @param k The column index.
 *
 * @return The 3D linear index with ik0 ordering.
 */
unsigned int
spamm_index_3D_ik0 (const unsigned int i, const unsigned int k)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 2;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if ((k & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;

    if ((i & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 2;
    getmask <<= 1;
  }

  return index;
}

/** Convert a 3D linear index to a 2D index.
 *
 * @param index_3D_i0j The 3D linear index.
 *
 * @return The 2D linear index.
 */
unsigned int
spamm_index_3D_i0j_to_2D (const unsigned int index_3D_i0j)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 1;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if ((index_3D_i0j & getmask) != 0)
    {
      index |= setmask;
    }
    getmask <<= 2;
    setmask <<= 1;

    if ((index_3D_i0j & getmask) != 0)
    {
      index |= setmask;
    }
    getmask <<= 1;
    setmask <<= 1;
  }

  return index;
}

/** Extract the k index from a 3D linear index.
 *
 * @param index_3D_ikj The 3D linear index in ikj encoding.
 *
 * @return The k index.
 */
unsigned int
spamm_index_3D_ikj_to_k (const unsigned int index_3D_ikj)
{
  short bit_index;
  unsigned int getmask = 2;
  unsigned int setmask = 1;
  unsigned int k = 0;

  for (bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if ((index_3D_ikj & getmask) != 0)
    {
      k |= setmask;
    }
    getmask <<= 3;
    setmask <<= 1;
  }

  return k;
}

/** Return a linear offset into a dense matrix block in row major order.
 *
 * @param i The row index.
 * @param j The column index.
 * @param M The number of rows of the matrix.
 * @param N The number of columns of the matrix.
 *
 * @return The linear offset into the matrix.
 */
unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N)
{
  assert(i < M);
  assert(j < N);

  return i*N+j;
}

/** Return a linear offset into a dense matrix block in column major order.
 *
 * @param i The row index.
 * @param j The column index.
 * @param M The number of rows of the matrix.
 * @param N The number of columns of the matrix.
 *
 * @return The linear offset into the matrix.
 */
unsigned int
spamm_index_column_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N)
{
  if (i >= M)
  {
    SPAMM_FATAL("i (%u) out of bounds (%u)\n", i, M);
  }

  if (j >= N)
  {
    SPAMM_FATAL("j (%u) out of bounds (%u)\n", j, N);
  }

  return i+j*M;
}

/** Return a linear offset in Z-curve ordering.
 *
 * @param i The row index.
 * @param j The column index.
 *
 * @return The offset.
 */
unsigned int
spamm_index_Z_curve (const unsigned int i, const unsigned int j)
{
  unsigned int i_internal = i;
  unsigned int j_internal = j;
  unsigned int setmask_i = 2;
  unsigned int setmask_j = 1;
  unsigned int offset = 0;

  while (i_internal > 0 || j_internal > 0)
  {
    if (i_internal > 0)
    {
      if ((i_internal & 0x1) != 0)
      {
        offset |= setmask_i;
      }
      i_internal >>= 1;
    }

    if (j_internal > 0)
    {
      if ((j_internal & 0x1) != 0)
      {
        offset |= setmask_j;
      }
      j_internal >>= 1;
    }

    setmask_i <<= 2;
    setmask_j <<= 2;
  }

  return offset;
}

/** Return a linear offset into the block norms at the kernel tier. The
 * indices are within the blocked kernel matrix, i.e. within the range of [0,
 * SPAMM_N_KERNEL_BLOCKED[.
 *
 * @param i The row index of the kernel block.
 * @param j The column index of the kernel block.
 *
 * @return The offset into the norm[] arry at the kernel tier.
 */
unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j)
{
  unsigned int offset;

  assert(i < SPAMM_N_KERNEL_BLOCKED);
  assert(j < SPAMM_N_KERNEL_BLOCKED);

  offset = spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCKED, SPAMM_N_KERNEL_BLOCKED);

  return offset;
}

/** Return a linear offset into a dense matrix block at the kernel tier. The
 * indices are within the kernel block at the kernel tier, i.e. within the
 * range of [0, SPAMM_N_KERNEL[.
 *
 * @param i The row index within the kernel block.
 * @param j The column index within the kernel block.
 * @param layout The layout of the basic matrix blocks at the kernel level.
 *
 * @return The offset into a kernel tier matrix block.
 */
unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout)
{
  assert(i < SPAMM_N_KERNEL);
  assert(j < SPAMM_N_KERNEL);

  return spamm_index_kernel_block_hierarchical(i/SPAMM_N_BLOCK, j/SPAMM_N_BLOCK, i%SPAMM_N_BLOCK, j%SPAMM_N_BLOCK, layout);
}

/** Return a linear offset into a dense matrix block at the kernel tier. The
 * indices are within the kernel block at the kernel tier, i.e. within the
 * range of [0, SPAMM_N_KERNEL[.
 *
 * @param i The row index within the kernel block.
 * @param j The column index within the kernel block.
 * @param layout The layout of the basic matrix blocks at the kernel level.
 *
 * @return The offset into a kernel tier matrix block.
 */
unsigned int
spamm_index_kernel_block_transpose (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout)
{
  assert(i < SPAMM_N_KERNEL);
  assert(j < SPAMM_N_KERNEL);

  return spamm_index_kernel_block_transpose_hierarchical(i/SPAMM_N_BLOCK, j/SPAMM_N_BLOCK, i%SPAMM_N_BLOCK, j%SPAMM_N_BLOCK, layout);
}

/** Return a linear offset into a kernel tier matrix block using hierarchical
 * indexing.
 *
 * @param i_blocked The row index of the basic matrix block, i.e. [0,
 * SPAMM_N_KERNEL_BLOCKED[.
 * @param j_blocked The column index of the basic matrix block, i.e. [0,
 * SPAMM_N_KERNEL_BLOCKED[.
 * @param i_basic The row index in the basic matrix block, i.e. [0,
 * SPAMM_N_BLOCK[.
 * @param j_basic The column index in the basic matrix block, i.e. [0,
 * SPAMM_N_BLOCK[.
 * @param layout The layout of the basic matrix blocks at the kernel level.
 *
 * @return The offset into the kernel matrix.
 */
unsigned int
spamm_index_kernel_block_hierarchical (const unsigned int i_blocked,
    const unsigned int j_blocked, const unsigned int i_basic,
    const unsigned int j_basic, const enum spamm_layout_t layout)
{
  unsigned int offset;

  assert(i_blocked < SPAMM_N_KERNEL_BLOCKED);
  assert(j_blocked < SPAMM_N_KERNEL_BLOCKED);
  assert(i_basic < SPAMM_N_BLOCK);
  assert(j_basic < SPAMM_N_BLOCK);

  switch (layout)
  {
    case dense_column_major:
      offset = spamm_index_column_major(i_blocked*SPAMM_N_BLOCK+i_basic, j_blocked*SPAMM_N_BLOCK+j_basic, SPAMM_N_KERNEL, SPAMM_N_KERNEL);
      break;

    case row_major:
      offset = spamm_index_row_major(i_basic, j_basic, SPAMM_N_BLOCK, SPAMM_N_BLOCK)
        +SPAMM_N_BLOCK*SPAMM_N_BLOCK * spamm_index_row_major(i_blocked, j_blocked, SPAMM_N_KERNEL_BLOCKED, SPAMM_N_KERNEL_BLOCKED);
      break;

    case Z_curve:
      offset = spamm_index_row_major(i_basic, j_basic, SPAMM_N_BLOCK, SPAMM_N_BLOCK)
        +SPAMM_N_BLOCK*SPAMM_N_BLOCK * spamm_index_Z_curve(i_blocked, j_blocked);
      break;

    default:
      SPAMM_FATAL("unknown layout (%i)\n", layout);
      break;
  }

  return offset;
}

/** Return a linear offset into a kernel tier matrix block using hierarchical
 * indexing. This function is used on the transpose part, i.e. the field
 * block_tranpose.
 *
 * @param i_blocked The row index of the basic matrix block, i.e. [0,
 * SPAMM_N_KERNEL_BLOCKED[.
 * @param j_blocked The column index of the basic matrix block, i.e. [0,
 * SPAMM_N_KERNEL_BLOCKED[.
 * @param i_basic The row index in the basic matrix block, i.e. [0,
 * SPAMM_N_BLOCK[.
 * @param j_basic The column index in the basic matrix block, i.e. [0,
 * SPAMM_N_BLOCK[.
 * @param layout The layout of the basic matrix blocks at the kernel level.
 *
 * @return The offset into the kernel matrix.
 */
unsigned int
spamm_index_kernel_block_transpose_hierarchical (const unsigned int i_blocked,
    const unsigned int j_blocked, const unsigned int i_basic,
    const unsigned int j_basic, const enum spamm_layout_t layout)
{
  unsigned int offset;

  assert(i_blocked < SPAMM_N_KERNEL_BLOCKED);
  assert(j_blocked < SPAMM_N_KERNEL_BLOCKED);
  assert(i_basic < SPAMM_N_BLOCK);
  assert(j_basic < SPAMM_N_BLOCK);

  switch (layout)
  {
    case dense_column_major:
      offset = spamm_index_column_major(j_blocked*SPAMM_N_BLOCK+j_basic, i_blocked*SPAMM_N_BLOCK+i_basic, SPAMM_N_KERNEL, SPAMM_N_KERNEL);
      break;

    case row_major:
      offset = spamm_index_row_major(j_basic, i_basic, SPAMM_N_BLOCK, SPAMM_N_BLOCK)
        +SPAMM_N_BLOCK*SPAMM_N_BLOCK * spamm_index_row_major(i_blocked, j_blocked, SPAMM_N_KERNEL_BLOCKED, SPAMM_N_KERNEL_BLOCKED);
      break;

    case Z_curve:
      offset = spamm_index_row_major(j_basic, i_basic, SPAMM_N_BLOCK, SPAMM_N_BLOCK)
        +SPAMM_N_BLOCK*SPAMM_N_BLOCK * spamm_index_Z_curve(i_blocked, j_blocked);
      break;

    default:
      SPAMM_FATAL("unknown layout (%i)\n", layout);
      break;
  }

  return offset;
}

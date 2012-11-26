#include "spamm.h"

#include <stdio.h>

#define BITFIELD_SIZE 15

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

  for(bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if((j & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;

    if((k & getmask) != 0)
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

  for(bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if((k & getmask) != 0)
    {
      index |= setmask;
    }
    setmask <<= 1;

    if((i & getmask) != 0)
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

  for(bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if((index_3D_i0j & getmask) != 0)
    {
      index |= setmask;
    }
    getmask <<= 2;
    setmask <<= 1;

    if((index_3D_i0j & getmask) != 0)
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

  for(bit_index = 0; bit_index < BITFIELD_SIZE; bit_index++)
  {
    if((index_3D_ikj & getmask) != 0)
    {
      k |= setmask;
    }
    getmask <<= 3;
    setmask <<= 1;
  }

  return k;
}

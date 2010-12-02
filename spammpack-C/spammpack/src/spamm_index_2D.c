#include "spamm.h"
#include <stdio.h>

#define BITFIELD_SIZE 10

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

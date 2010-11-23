#include "spamm.h"
#include <stdio.h>

#define BITFIELD_SIZE 15

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

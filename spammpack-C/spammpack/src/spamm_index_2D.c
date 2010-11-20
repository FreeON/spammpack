#include "spamm.h"

#include <stdio.h>

unsigned int
spamm_index_2D (const unsigned int i, const unsigned int j)
{
  short bit_index;
  unsigned int getmask = 1;
  unsigned int setmask = 1;
  unsigned int index = 0;

  for (bit_index = 0; bit_index < 10; bit_index++)
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

  printf("(%u,%u) --> %u\n", i, j, index);
  return index;
}

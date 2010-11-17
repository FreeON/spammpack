#include <stdio.h>

/* Generate the 2-D dilated integer table from
 *
 * IEEE Trans. Comput. C99, 2, R. Raman and D. S. Wise (2007).
 */
int
main ()
{
  unsigned short bitmask, result_bitmask;
  unsigned short x, y, i, result;

  for (x = 0; x <= 255; x += 8)
  {
    for (y = x; y < x+8; y++)
    {
      result = 0;
      bitmask = 1;
      result_bitmask = 1;

      for (i = 0; i < 8; i++)
      {
        if (x & bitmask)
        {
          result |= result_bitmask;
        }
        bitmask <<= 1;
        result_bitmask <<= 2;
      }
      printf("0x%04x, ", result);
    }
    printf("\n");
  }

  return 0;
}

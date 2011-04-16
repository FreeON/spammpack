#include "spamm.h"
#include <assert.h>

/** Convert an unsigned int to a binary string.
 *
 * @param width The number of bits to convert.
 * @param i The integer to convert.
 * @param result The resulting string. This string has to be allocated by the
 * caller and has to be big enough. This function does not check whether
 * that's true.
 */
void
spamm_uint_to_bin_string (const unsigned int width, const unsigned int i, char *result)
{
  short index;
  unsigned int bitmask = 1 << (width-1);

  assert(result != NULL);

  if (width == 0)
  {
    result[0] = '0';
    result[1] = '\0';
  }

  else
  {
    for (index = width-1; index >= 0; index--)
    {
      if ((i & bitmask) == 0)
      {
        result[width-1-index] = '0';
      }

      else
      {
        result[width-1-index] = '1';
      }

      bitmask >>= 1;
    }
    result[width-1-index] = '\0';
  }
}

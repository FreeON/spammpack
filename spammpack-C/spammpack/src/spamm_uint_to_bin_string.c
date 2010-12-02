#include "spamm.h"
#include <assert.h>

#define STRING_WIDTH 10

/** Convert an unsigned int to a binary string.
 *
 * @param i The integer to convert.
 * @param result The resulting string. This string has to be allocated by the
 * caller and has to be big enough. This function does not check whether
 * that's true.
 */
void
spamm_uint_to_bin_string (const unsigned int i, char *result)
{
  short index;
  unsigned int bitmask = 1 << (STRING_WIDTH-1);

  assert(result != NULL);

  for (index = STRING_WIDTH-1; index >= 0; index--)
  {
    if ((i & bitmask) == 0)
    {
      result[STRING_WIDTH-1-index] = '0';
    }

    else
    {
      result[STRING_WIDTH-1-index] = '1';
    }

    bitmask >>= 1;
  }
  result[STRING_WIDTH-1-index] = '\0';
}

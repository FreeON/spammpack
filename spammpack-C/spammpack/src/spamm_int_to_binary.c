#include "spamm.h"

/** Convert an integer to string in binary.
 *
 * The parameter binary_string has to be pre-allocated by the caller and allow
 * for enough room to fit the output string. This function does not test for
 * whether that is true.
 *
 * @param integer The integer to convert.
 * @param width The width of the string, i.e. how many bits should be printed.
 * @param binary_string The string.
 */
void
spamm_int_to_binary (const unsigned int integer, const int width, char *binary_string)
{
  int i;
  unsigned int mask = 1;

  if (width == 0)
  {
    binary_string[0] = '0';
    binary_string[1] = '\0';
  }

  else
  {
    for (i = 0; i < width; ++i)
    {
      binary_string[i] = '0';
    }
    binary_string[width] = '\0';

    for (i = 0; i < width; ++i)
    {
      if (mask & integer) { binary_string[width-1-i] = '1'; }
      mask = mask << 1;
    }
  }
}

#include "spamm.h"

/** Implementation of the fast exponention algorithm by squaring. This
 * function calculates for integer arguments \f$ b^n /f$.
 *
 * @param b The base \f$ b \f$.
 * @param n The exponent \f$ n \f$
 *
 * @return \f$ b^n /f$.
 */
unsigned int
ipow (unsigned int b, unsigned int n)
{
  int result = 1;

  if (b == 2)
  {
    result <<= n;
  }

  else
  {
    while (n != 0)
    {
      if ((n & 1) != 0) /* Odd exponent, multiply by b. */
      {
        result *= b;
      }
      n >>= 1; /* Divide n by 2. */
      b *= b;  /* Square b. */
    }
  }

  return result;
}

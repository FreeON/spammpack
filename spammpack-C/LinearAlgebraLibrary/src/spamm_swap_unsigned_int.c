#include "spamm.h"

/** Swap two unsigned int.
 *
 * @param a The first unsigned int.
 * @param b The second unsigned int.
 */
void
spamm_swap_unsigned_int (unsigned int *a, unsigned int *b)
{
  unsigned int temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

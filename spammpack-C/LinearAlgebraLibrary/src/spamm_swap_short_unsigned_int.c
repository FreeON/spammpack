#include "spamm.h"

/** Swap two short unsigned int.
 *
 * @param a The first short unsigned int.
 * @param b The second short unsigned int.
 */
void
spamm_swap_short_unsigned_int (short unsigned int *a, short unsigned int *b)
{
  short unsigned int temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

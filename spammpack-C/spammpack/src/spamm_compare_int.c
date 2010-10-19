#include "spamm.h"

/** Compare to integers.
 *
 * @param integer1 The first integer.
 * @param integer2 The second integer.
 *
 * @return integer1 < integer2: -1, integer1 == integer2: 0, integer1 >
 *         integer2: +1.
 */
int
spamm_compare_int (const void *integer1, const void *integer2)
{
  const int *int1 = integer1;
  const int *int2 = integer2;

  if (*int1 < *int2) { return -1; }
  else if (*int1 == *int2) { return 0; }
  else { return 1; }
}

#include "spamm.h"

/** Swap to floating_point_t.
 *
 * @param a The first floating_point_t.
 * @param b The second floating_point_t.
 */
void
spamm_swap_floating_point_t (floating_point_t *a, floating_point_t *b)
{
  floating_point_t temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

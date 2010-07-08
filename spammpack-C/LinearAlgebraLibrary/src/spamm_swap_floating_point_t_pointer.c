#include "spamm.h"

/** Swap two pointers to a floating_point_t.
 *
 * @param a A pointer to the first floating_point_t.
 * @param b A pointer to the second floating_point_t.
 */
void
spamm_swap_floating_point_t_pointer (floating_point_t **a, floating_point_t **b)
{
  floating_point_t *temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

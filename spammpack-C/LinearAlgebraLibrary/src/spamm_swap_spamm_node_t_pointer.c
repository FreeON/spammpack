#include "spamm.h"

/** Swap two pointers to a spamm_node_t.
 *
 * @param a A pointer to the first spamm_node_t.
 * @param b A pointer to the second spamm_node_t.
 */
void
spamm_swap_spamm_node_t_pointer (struct spamm_node_t **a, struct spamm_node_t **b)
{
  struct spamm_node_t *temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

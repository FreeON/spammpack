#include "spamm.h"

/** Compare to unsigned int values. This is a helper function for the
 * hashtables in dealing with the linear indices.
 *
 * @param a A pointer to the first index.
 * @param b A pointer to the second index.
 *
 * @return TRUE if the two values are equal, FALSE if not.
 */
gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b)
{
  const unsigned int *value_1 = a;
  const unsigned int *value_2 = b;

  if ((*value_1) == (*value_2)) { return TRUE; }
  else { return FALSE; }
}

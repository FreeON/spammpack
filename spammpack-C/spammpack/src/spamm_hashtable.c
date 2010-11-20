#include "spamm.h"

gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b)
{
  const unsigned int *value_1 = a;
  const unsigned int *value_2 = b;

  if ((*value_1) == (*value_2)) { return TRUE; }
  else { return FALSE; }
}
